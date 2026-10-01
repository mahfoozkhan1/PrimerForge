import os
import time
import xml.etree.ElementTree as ET

import requests


BASE = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def _headers(tool: str, email: str):
    return {
        "User-Agent": f"{tool} {email}".strip(),
        "Accept": "text/plain",
    }


def _extract_value(text: str, key: str):
    """
    NCBI BLAST can return values in slightly different response formats.
    Try several common formats instead of assuming only `RID = ...`.
    """
    for line in text.splitlines():
        line = line.strip()

        if line.startswith(f"{key} ="):
            return line.split("=", 1)[1].strip()

        if line.startswith(f"{key}="):
            return line.split("=", 1)[1].strip()

        if line.startswith(f"{key}\t"):
            return line.split("\t", 1)[1].strip()

    return None


def _extract_rid(text: str):
    rid = _extract_value(text, "RID")

    if rid:
        return rid

    # Some responses may contain the RID in a less structured form.
    for line in text.splitlines():
        if "RID" in line.upper():
            parts = line.replace(":", "=").split("=", 1)
            if len(parts) == 2:
                candidate = parts[1].strip()

                if candidate and len(candidate) >= 5:
                    return candidate

    return None


def _submit_blast(
    query: str,
    db: str,
    organism: str | None,
    email: str,
    tool: str,
    api_key: str = "",
):
    """
    Submit one BLASTN job to NCBI.

    Primer sequences are short, so SHORT_QUERY_ADJUST is enabled.
    """

    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": db,
        "QUERY": query,
        "MEGABLAST": "on",
        "SHORT_QUERY_ADJUST": "on",
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": "20",
        "TOOL": tool,
        "EMAIL": email,
    }

    if organism:
        params["ENTREZ_QUERY"] = f'"{organism}"[Organism]'

    if api_key:
        params["API_KEY"] = api_key

    response = requests.post(
        BASE,
        data=params,
        headers=_headers(tool, email),
        timeout=60,
    )

    response.raise_for_status()

    body = response.text

    rid = _extract_rid(body)

    if not rid:
        # Preserve NCBI's actual response so the UI can show the
        # real reason instead of only "RID not returned".
        preview = body.strip()

        if len(preview) > 1000:
            preview = preview[:1000] + "..."

        raise RuntimeError(
            "NCBI BLAST did not return a request ID. "
            f"NCBI response: {preview}"
        )

    return rid


def _wait_for_blast(
    rid: str,
    email: str,
    tool: str,
    timeout: int = 180,
):
    """
    Poll NCBI until the BLAST job is ready.

    We deliberately use a conservative polling interval because
    NCBI BLAST is a shared service.
    """

    start = time.time()

    while time.time() - start < timeout:

        time.sleep(10)

        response = requests.get(
            BASE,
            params={
                "CMD": "Get",
                "RID": rid,
                "FORMAT_OBJECT": "SearchInfo",
            },
            headers=_headers(tool, email),
            timeout=60,
        )

        response.raise_for_status()

        body = response.text

        if "Status=READY" in body:
            return

        if "Status=FAILED" in body:
            raise RuntimeError(
                f"NCBI BLAST job failed for RID {rid}."
            )

        if "Status=UNKNOWN" in body:
            raise RuntimeError(
                f"NCBI BLAST no longer recognizes RID {rid}."
            )

    raise TimeoutError(
        f"NCBI BLAST timed out after {timeout} seconds."
    )


def _get_results(
    rid: str,
    email: str,
    tool: str,
):
    response = requests.get(
        BASE,
        params={
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "XML",
        },
        headers=_headers(tool, email),
        timeout=90,
    )

    response.raise_for_status()

    return summarize_xml(response.text)


def run_blast(
    query: str,
    db="nt",
    organism=None,
    timeout=180,
    email=None,
    tool=None,
    api_key=None,
):
    """
    Run a complete NCBI BLASTN search.

    This function is intentionally conservative with NCBI requests.
    """

    email = email or os.getenv("NCBI_EMAIL", "").strip()

    tool = tool or os.getenv(
        "NCBI_TOOL",
        "PrimerForge",
    ).strip()

    api_key = api_key or os.getenv(
        "NCBI_API_KEY",
        "",
    ).strip()

    if not email:
        raise RuntimeError(
            "NCBI_EMAIL is not configured. "
            "Add NCBI_EMAIL to the Render environment variables."
        )

    query = "".join(
        query.split()
    ).upper()

    if len(query) < 10:
        raise ValueError(
            "BLAST query is too short for reliable primer specificity analysis."
        )

    rid = _submit_blast(
        query=query,
        db=db,
        organism=organism,
        email=email,
        tool=tool,
        api_key=api_key,
    )

    _wait_for_blast(
        rid=rid,
        email=email,
        tool=tool,
        timeout=timeout,
    )

    return _get_results(
        rid=rid,
        email=email,
        tool=tool,
    )


def summarize_xml(xml: str):
    """
    Convert the BLAST XML result into a compact structure
    suitable for PrimerForge ranking and display.
    """

    root = ET.fromstring(xml)

    hits = []

    for hit in root.findall(".//Hit")[:20]:

        hsps = hit.findall("./Hit_hsps/Hsp")

        if not hsps:
            continue

        hsp = hsps[0]

        identity = int(
            hsp.findtext(
                "Hsp_identity",
                "0",
            )
        )

        alignment_length = int(
            hsp.findtext(
                "Hsp_align-len",
                "1",
            )
        )

        gaps = int(
            hsp.findtext(
                "Hsp_gaps",
                "0",
            )
        )

        hits.append(
            {
                "accession": hit.findtext(
                    "Hit_accession"
                ),
                "title": hit.findtext(
                    "Hit_def"
                ),
                "identity_percent": round(
                    100 * identity / alignment_length,
                    2,
                ),
                "alignment_length": alignment_length,
                "gaps": gaps,
                "evalue": float(
                    hsp.findtext(
                        "Hsp_evalue",
                        "0",
                    )
                ),
                "bit_score": float(
                    hsp.findtext(
                        "Hsp_bit-score",
                        "0",
                    )
                ),
            }
        )

    return hits
