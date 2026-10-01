# PrimerForge — Automated Primer Design & Validation

PrimerForge is a modular web application for automated oligonucleotide design and validation. It accepts a nucleotide FASTA sequence, lets the user select an assay type, designs candidates with Primer3, performs thermodynamic checks (Tm, GC%, hairpin, homodimer/self-complementarity and heterodimer/pair complementarity), optionally runs NCBI BLAST specificity checks, ranks candidates, and produces a machine-readable JSON report.

> **Important:** This is a research/engineering tool, not a substitute for experimental validation. Primer design rules depend on polymerase, chemistry, sample type, organism, assay platform, and laboratory SOPs. BLAST specificity depends strongly on the selected database and organism filter.

## Supported modes

- Conventional PCR
- RT-PCR / cDNA amplification
- qPCR
- ddPCR
- miRNA assay (sequence-specific forward oligo + universal/stem-loop RT workflow guidance)
- siRNA candidate design (21-nt duplex candidates; off-target BLAST is optional)

## What it checks

For DNA primer candidates:

- Length
- GC percentage
- Tm
- 3' end GC/clamp heuristic
- Homopolymer runs
- Hairpin thermodynamics
- Homodimer thermodynamics
- Heterodimer thermodynamics
- Primer3 self-complementarity / 3' self-complementarity
- Primer3 pair complementarity
- Amplicon size
- Optional NCBI BLAST specificity

For siRNA candidates:

- 21-nt guide candidates
- GC percentage
- Homopolymer runs
- 3'/5' composition heuristics
- Simple seed-risk flagging
- Optional BLAST screening

## Architecture

```text
Browser
  |
  v
FastAPI (/api/analyze)
  |
  +--> Input validation / FASTA parsing
  +--> Assay-specific parameter preset
  +--> Primer3 candidate generation
  +--> Thermodynamic validation
  +--> Optional NCBI BLAST specificity
  +--> Ranking engine
  +--> JSON report
```

## Local setup

### 1. Python

Python 3.10+ recommended.

### 2. Install dependencies

```bash
python -m venv .venv
# Windows: .venv\\Scripts\\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

### 3. Start

```bash
uvicorn app.main:app --reload
```

Open `http://127.0.0.1:8000`.

## BLAST

The application uses NCBI's public BLAST URL API when BLAST is enabled. NCBI provides an HTTPS BLAST API for submitting jobs and retrieving results. Keep request volume conservative and identify the application with a descriptive User-Agent/email in `.env`. For higher-volume institutional workflows, use a local BLAST+ installation and local database instead.

Create `.env` from `.env.example`:

```env
NCBI_EMAIL=you@example.com
NCBI_TOOL=PrimerForge/1.0
NCBI_API_KEY=
BLAST_DB=nt
```

The UI lets you choose whether to run BLAST. For a human qPCR assay, select `Homo sapiens` and a suitable nucleotide database. For a different organism, use the organism field.

## GitHub Actions

The included workflow runs tests on every push/PR. The application itself is designed to run on a free Python host or locally; NCBI BLAST calls require outbound internet access.

## API

`POST /api/analyze`

Example payload:

```json
{
  "assay_type": "qpcr",
  "sequence": ">TP53\nATGG...",
  "organism": "Homo sapiens",
  "run_blast": false,
  "max_pairs": 10
}
```

## Scientific implementation notes

Primer3 is used for conventional primer-pair generation and provides thermodynamic calculations for Tm, hairpins, homodimers and heterodimers. Primer3's own manual documents the self-complementarity, 3' self-complementarity, hairpin and pair-complementarity metrics used by this application.

The default presets are intentionally conservative starting points, not universal laboratory standards. Edit `app/presets.py` for your exact chemistry/SOP.


---

**Made by Md Mahfooz Khan**

PrimerForge is an in-silico design and validation aid; experimental confirmation and assay-specific protocol review remain essential.
