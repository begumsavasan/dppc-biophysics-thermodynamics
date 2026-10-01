# DPPC Biophysics & Thermodynamics

Public companion repository for reproducible work on hydrated DPPC membrane phase behavior, thermodynamic structure, and pressure–temperature analysis.

## Archived research output

**Mendeley Data, Version 1**  
DOI: https://doi.org/10.17632/r97b37bfjg.1

The archived package was published on 17 September 2026 and contains research-case records, evidence mapping, provenance-labelled transcribed source-table inputs, dated reconstruction code and outputs, a preserved historical implementation/test suite, and a DPPC source audit.

The archive reports **no new experimental measurements**. Its case records should not be interpreted as independent experiments or independent validations.

## Scientific scope

This work addresses:

- hydrated-DPPC membrane phase behavior
- thermodynamic consistency and pressure–temperature relationships
- phase-boundary geometry
- equation-based and numerical analysis
- estimator and null-model auditing
- provenance-aware scientific computing
- reproducible reconstruction of archived calculations

The associated research record does not claim that the hydrated-DPPC Lβ′/LβI boundary has been assigned a termination coordinate. The public archive is structured to preserve that distinction rather than overstate a result.

## Executable reconstruction subset

For lightweight inspection without downloading the complete archive, this repository mirrors a small executable subset from the canonical Mendeley Data package:

- [reconstruction/data/source_rows.csv](reconstruction/data/source_rows.csv)
- [reconstruction/run_analysis.py](reconstruction/run_analysis.py)
- [reconstruction/verify_reconstruction.py](reconstruction/verify_reconstruction.py)
- [reconstruction/README.md](reconstruction/README.md)

The mirrored reconstruction regenerates the archived analysis outputs and the verification script completes **47 numerical regression checks**. These are regression checks on reconstructed calculations, not 47 independent scientific validations.

## Reproduce the complete archived package

The complete reproducibility package remains distributed through Mendeley Data.

Recommended environment:

- Python 3.12
- NumPy 2.3.5
- SciPy 1.17.0
- pytest 9.1.1

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the full public reproduction map.

## Repository boundary

This repository is a **public-facing companion**, not the private manuscript/reviewer repository.

It does not contain:

- unpublished manuscript text
- reviewer correspondence
- confidential submission material
- restricted source-publication PDFs
- private collaboration records
- unreleased intellectual property

Public material is limited to documentation and content already suitable for open release.

## Citation

See [CITATION.cff](CITATION.cff).

## Researcher

Begüm Savaşan Asgarlı  
ORCID: https://orcid.org/0000-0003-4809-7115
