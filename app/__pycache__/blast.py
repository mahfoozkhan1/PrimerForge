import os
import time
import requests
from urllib.parse import urlencode

BASE = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def run_blast(query: str, db="nt", organism=None, timeout=120, email=None, tool=None, api_key=None):
    email = email or os.getenv("NCBI_EMAIL", "")
    tool = tool or os.getenv("NCBI_TOOL", "PrimerForge/1.0")
    api_key = api_key or os.getenv("NCBI_API_KEY", "")
    headers = {"User-Agent": f"{tool} {email}".strip()}
    params = {
        "CMD": "Put", "PROGRAM": "blastn", "DATABASE": db, "QUERY": query,
        "MEGABLAST": "on", "FORMAT_TYPE": "XML", "HITLIST_SIZE": "20",
        "TOOL": tool, "EMAIL": email,
    }
    if organism:
        params["ENTREZ_QUERY"] = f'"{organism}"[Organism]'
    if api_key:
        params["API_KEY"] = api_key
    r = requests.post(BASE, data=params, headers=headers, timeout=30)
    r.raise_for_status()
    text = r.text
    rid = _extract(text, "RID")
    if not rid:
        raise RuntimeError("NCBI BLAST did not return a request ID.")
    start = time.time()
    while time.time() - start < timeout:
        time.sleep(4)
        status = requests.get(BASE, params={"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"}, headers=headers, timeout=30)
        status.raise_for_status()
        if "Status=READY" in status.text:
            break
        if "Status=FAILED" in status.text:
            raise RuntimeError("NCBI BLAST job failed.")
    else:
        raise TimeoutError("NCBI BLAST timed out.")
    result = requests.get(BASE, params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}, headers=headers, timeout=60)
    result.raise_for_status()
    return summarize_xml(result.text)


def _extract(text, key):
    for line in text.splitlines():
        if line.startswith(key + " ="):
            return line.split("=", 1)[1].strip()
    return None


def summarize_xml(xml: str):
    # Lightweight parser: enough for specificity ranking without storing the entire XML.
    import xml.etree.ElementTree as ET
    root = ET.fromstring(xml)
    hits = []
    for hit in root.findall(".//Hit")[:20]:
        hsps = hit.findall("./Hit_hsps/Hsp")
        if not hsps:
            continue
        h = hsps[0]
        ident = int(h.findtext("Hsp_identity", "0"))
        align = int(h.findtext("Hsp_align-len", "1"))
        gaps = int(h.findtext("Hsp_gaps", "0"))
        hits.append({
            "accession": hit.findtext("Hit_accession"),
            "title": hit.findtext("Hit_def"),
            "identity_percent": round(100 * ident / align, 2),
            "alignment_length": align,
            "gaps": gaps,
            "evalue": float(h.findtext("Hsp_evalue", "0")),
            "bit_score": float(h.findtext("Hsp_bit-score", "0")),
        })
    return hits
