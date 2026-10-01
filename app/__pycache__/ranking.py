
def score_pair(pair, blast_forward=None, blast_reverse=None):
    score = 100.0
    reasons = []
    # Balanced Tm and GC.
    for label, tm in [("F", pair["forward_tm"]), ("R", pair["reverse_tm"])]:
        score -= min(15, abs(float(tm) - 60.0) * 1.8)
    score -= min(10, abs(float(pair["forward_gc"]) - 50.0) * 0.35)
    score -= min(10, abs(float(pair["reverse_gc"]) - 50.0) * 0.35)
    score -= min(8, abs(float(pair["forward_thermo"]["heterodimer"]["tm"]) if "heterodimer" in pair["forward_thermo"] else 0)) if False else 0
    hd_tm = float(pair["heterodimer"]["tm"])
    if pair["heterodimer"]["structure_found"] and hd_tm >= 45:
        score -= 15
        reasons.append("strong heterodimer signal")
    for side in ("forward_thermo", "reverse_thermo"):
        hp = pair[side]["hairpin"]
        if hp["structure_found"] and hp["tm"] >= 45:
            score -= 10
            reasons.append("hairpin risk")
    for side in ("forward_thermo", "reverse_thermo"):
        hd = pair[side]["homodimer"]
        if hd["structure_found"] and hd["tm"] >= 45:
            score -= 8
            reasons.append("homodimer risk")
    for key in ("forward_homopolymer_max", "reverse_homopolymer_max"):
        if pair["heuristics"][key] >= 5:
            score -= 8
            reasons.append("homopolymer risk")
    if blast_forward is not None:
        # Penalize if multiple strong hits exist; exact interpretation depends on assay/organism.
        strong = sum(h["identity_percent"] >= 90 and h["alignment_length"] >= 18 for h in blast_forward)
        if strong > 1:
            score -= min(20, (strong - 1) * 5)
            reasons.append("multiple high-identity BLAST hits for forward primer")
    if blast_reverse is not None:
        strong = sum(h["identity_percent"] >= 90 and h["alignment_length"] >= 18 for h in blast_reverse)
        if strong > 1:
            score -= min(20, (strong - 1) * 5)
            reasons.append("multiple high-identity BLAST hits for reverse primer")
    pair["quality_score"] = round(max(0.0, score), 2)
    pair["ranking_notes"] = reasons or ["No major automated red flags"]
    return pair
