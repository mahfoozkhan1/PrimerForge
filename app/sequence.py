import re
from Bio.Seq import Seq

DNA_RE = re.compile(r"^[ACGTUNRYKMSWBDHV]+$", re.I)

def parse_sequence(raw: str):
    lines = [x.strip() for x in raw.strip().splitlines() if x.strip()]
    header = "sequence"
    if lines and lines[0].startswith(">"):
        header = lines[0][1:].strip() or "sequence"
        lines = lines[1:]
    seq = "".join(lines).upper().replace("U", "T").replace(" ", "")
    if not seq:
        raise ValueError("No nucleotide sequence was supplied.")
    if not DNA_RE.match(seq):
        bad = sorted(set(re.sub(r"[ACGTNRYKMSWBDHV]", "", seq)))
        raise ValueError(f"Invalid nucleotide characters: {bad}")
    return header, seq

def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def gc_percent(seq: str) -> float:
    informative = [b for b in seq if b in "ACGT"]
    return round(100.0 * sum(b in "GC" for b in informative) / len(informative), 2) if informative else 0.0

def homopolymer_max(seq: str) -> int:
    best = cur = 1
    for a, b in zip(seq, seq[1:]):
        cur = cur + 1 if a == b else 1
        best = max(best, cur)
    return best
