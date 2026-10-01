import os
from pathlib import Path
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from dotenv import load_dotenv
from .models import AnalyzeRequest
from .sequence import parse_sequence, gc_percent
from .designer import design_pairs, design_sirna
from .ranking import score_pair
from .blast import run_blast

load_dotenv()
app = FastAPI(title="PrimerForge", version="1.0.0")
app.mount("/static", StaticFiles(directory=Path(__file__).parent.parent / "static"), name="static")

@app.get("/", response_class=HTMLResponse)
def home():
    return (Path(__file__).parent.parent / "static" / "index.html").read_text(encoding="utf-8")

@app.get("/api/health")
def health():
    return {"status": "ok", "service": "PrimerForge"}

@app.post("/api/analyze")
def analyze(req: AnalyzeRequest):
    try:
        header, seq = parse_sequence(req.sequence)
    except ValueError as e:
        raise HTTPException(400, str(e))
    min_len = 21 if req.assay_type == "sirna" else (18 if req.assay_type == "mirna" else 40)
    if len(seq) < min_len:
        raise HTTPException(400, f"Sequence is too short for {req.assay_type}; minimum input length is {min_len} nt.")

    warnings = []
    if req.assay_type == "sirna":
        candidates = design_sirna(seq, req.max_pairs * 2)
        if req.run_blast:
            for c in candidates[:req.max_pairs]:
                try:
                    c["blast"] = run_blast(c["sense_5to3"], req.blast_db, req.organism)
                except Exception as e:
                    warnings.append(f"BLAST failed for siRNA candidate {c['position']}: {e}")
        best = candidates[0] if candidates else None
        return {"status": "completed", "assay_type": req.assay_type, "sequence_id": header, "sequence_length": len(seq), "input_gc_percent": gc_percent(seq), "candidates": candidates, "best": best, "warnings": warnings}

    pairs, explain = design_pairs(seq, req.assay_type if req.assay_type != "mirna" else "mirna", req.max_pairs)
    if not pairs:
        warnings.append("Primer3 returned no primer pairs under the current constraints.")
    for pair in pairs:
        if req.run_blast:
            try:
                pair["blast_forward"] = run_blast(pair["forward"], req.blast_db, req.organism)
                pair["blast_reverse"] = run_blast(pair["reverse"], req.blast_db, req.organism)
            except Exception as e:
                warnings.append(f"BLAST failed for pair {pair['rank_from_primer3']}: {e}")
        score_pair(pair, pair.get("blast_forward"), pair.get("blast_reverse"))
    pairs.sort(key=lambda x: (-x.get("quality_score", 0), x.get("primer3_penalty", 9999)))
    best = pairs[0] if pairs else None
    if req.assay_type == "mirna":
        warnings.append("miRNA workflows are assay-specific. The pair shown is a sequence-specific DNA oligo design; validated stem-loop RT or poly(A)-tailing chemistry should be selected according to the kit/protocol used.")
    return {"status": "completed", "assay_type": req.assay_type, "sequence_id": header, "sequence_length": len(seq), "input_gc_percent": gc_percent(seq), "primer3_explain": explain, "candidates": pairs, "best": best, "warnings": warnings}
