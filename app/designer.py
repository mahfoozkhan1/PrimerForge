import primer3
from .presets import PRESETS
from .sequence import gc_percent, homopolymer_max, revcomp
from .thermo import thermo_metrics, heterodimer_metrics


def _primer3_params(p):
    return {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_NUM_RETURN": 30,
        "PRIMER_MIN_SIZE": p["primer_min_size"],
        "PRIMER_OPT_SIZE": p["primer_opt_size"],
        "PRIMER_MAX_SIZE": p["primer_max_size"],
        "PRIMER_MIN_TM": p["primer_min_tm"],
        "PRIMER_OPT_TM": p["primer_opt_tm"],
        "PRIMER_MAX_TM": p["primer_max_tm"],
        "PRIMER_MIN_GC": p["primer_min_gc"],
        "PRIMER_MAX_GC": p["primer_max_gc"],
        "PRIMER_PRODUCT_SIZE_RANGE": p["product_size_range"],
        "PRIMER_MAX_SELF_ANY": p["max_self_any"],
        "PRIMER_MAX_SELF_END": p["max_self_end"],
        "PRIMER_PAIR_MAX_COMPL_ANY": p["max_pair_any"],
        "PRIMER_PAIR_MAX_COMPL_END": p["max_pair_end"],
        "PRIMER_MAX_HAIRPIN_TH": p["max_hairpin_tm"],
        "PRIMER_EXPLAIN_FLAG": 1,
    }


def design_pairs(seq: str, assay_type: str, max_pairs: int = 10):
    p = PRESETS[assay_type]
    result = primer3.bindings.design_primers(
        {"SEQUENCE_ID": "target", "SEQUENCE_TEMPLATE": seq},
        _primer3_params(p),
    )
    pairs = []
    for i in range(min(max_pairs, int(result.get("PRIMER_PAIR_NUM_RETURNED", 0)))):
        l = result[f"PRIMER_LEFT_{i}_SEQUENCE"]
        r = result[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        pair = {
            "rank_from_primer3": i + 1,
            "forward": l,
            "reverse": r,
            "forward_start": result.get(f"PRIMER_LEFT_{i}", [None, None])[0],
            "reverse_start": result.get(f"PRIMER_RIGHT_{i}", [None, None])[0],
            "amplicon_size": result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
            "primer3_penalty": result.get(f"PRIMER_PAIR_{i}_PENALTY"),
            "forward_tm": result.get(f"PRIMER_LEFT_{i}_TM"),
            "reverse_tm": result.get(f"PRIMER_RIGHT_{i}_TM"),
            "forward_gc": result.get(f"PRIMER_LEFT_{i}_GC_PERCENT"),
            "reverse_gc": result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT"),
            "forward_self_any": result.get(f"PRIMER_LEFT_{i}_SELF_ANY"),
            "reverse_self_any": result.get(f"PRIMER_RIGHT_{i}_SELF_ANY"),
            "forward_self_end": result.get(f"PRIMER_LEFT_{i}_SELF_END"),
            "reverse_self_end": result.get(f"PRIMER_RIGHT_{i}_SELF_END"),
            "pair_compl_any": result.get(f"PRIMER_PAIR_{i}_COMPL_ANY"),
            "pair_compl_end": result.get(f"PRIMER_PAIR_{i}_COMPL_END"),
            "forward_hairpin_tm": result.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH"),
            "reverse_hairpin_tm": result.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH"),
        }
        pair["forward_thermo"] = thermo_metrics(l)
        pair["reverse_thermo"] = thermo_metrics(r)
        pair["heterodimer"] = heterodimer_metrics(l, r)
        pair["heuristics"] = {
            "forward_homopolymer_max": homopolymer_max(l),
            "reverse_homopolymer_max": homopolymer_max(r),
            "forward_3prime_base": l[-1],
            "reverse_3prime_base": r[-1],
        }
        pairs.append(pair)
    return pairs, result.get("PRIMER_PAIR_EXPLAIN")


def design_sirna(seq: str, max_candidates=20):
    candidates = []
    for i in range(0, len(seq) - 20):
        sense = seq[i:i+21]
        gc = gc_percent(sense)
        if 30 <= gc <= 55 and homopolymer_max(sense) <= 4:
            antisense = revcomp(sense)
            # A transparent heuristic, not a replacement for an established siRNA scoring model.
            score = 100 - abs(gc - 42.5) * 1.5 - (homopolymer_max(sense) - 2) * 4
            candidates.append({"position": i + 1, "sense_5to3": sense, "antisense_5to3": antisense, "gc_percent": gc, "score": round(score, 2), "seed_region": sense[1:8]})
    candidates.sort(key=lambda x: x["score"], reverse=True)
    return candidates[:max_candidates]
