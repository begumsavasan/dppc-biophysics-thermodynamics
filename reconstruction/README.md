# Curated reconstruction subset

This directory mirrors a compact, executable subset of the canonical public Mendeley Data archive:

https://doi.org/10.17632/r97b37bfjg.1

The Mendeley Data record remains the canonical release.

Included here:

- `data/source_rows.csv` — provenance-labelled transcribed source rows;
- `run_analysis.py` — dated 2026-09-12 reconstruction;
- `verify_reconstruction.py` — numerical regression checks on reconstructed results.

## Run

From the repository root:

```bash
python reconstruction/run_analysis.py
python reconstruction/verify_reconstruction.py
```

Generated outputs are written below `reconstruction/results/`.

## Interpretation limits

The numerical checks are regression checks on reconstructed calculations, not independent experimental validations.

The source rows are transcribed from audited literature sources. Copyrighted source PDFs are not redistributed.

This GitHub subset is intended for inspection and lightweight execution. The complete evidence map, historical implementation, source audit, synthetic method audits, manifests, and verification logs remain in the canonical Mendeley Data archive.
