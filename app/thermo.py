import primer3


def thermo_metrics(seq: str, mv=50.0, dv=1.5, dntp=0.6, dna=50.0):
    tm = primer3.calc_tm(seq, mv_conc=mv, dv_conc=dv, dntp_conc=dntp, dna_conc=dna)
    hp = primer3.calc_hairpin(seq, mv_conc=mv, dv_conc=dv, dntp_conc=dntp, dna_conc=dna)
    hd = primer3.calc_homodimer(seq, mv_conc=mv, dv_conc=dv, dntp_conc=dntp, dna_conc=dna)
    return {
        "tm": round(float(tm), 2),
        "hairpin": {"structure_found": bool(hp.structure_found), "tm": round(float(hp.tm), 2), "dg": round(float(hp.dg), 2)},
        "homodimer": {"structure_found": bool(hd.structure_found), "tm": round(float(hd.tm), 2), "dg": round(float(hd.dg), 2)},
    }


def heterodimer_metrics(a: str, b: str, mv=50.0, dv=1.5, dntp=0.6, dna=50.0):
    hd = primer3.calc_heterodimer(a, b, mv_conc=mv, dv_conc=dv, dntp_conc=dntp, dna_conc=dna)
    return {"structure_found": bool(hd.structure_found), "tm": round(float(hd.tm), 2), "dg": round(float(hd.dg), 2)}
