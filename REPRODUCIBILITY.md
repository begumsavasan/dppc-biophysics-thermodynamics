# Reproducibility

The canonical public reproducibility package for this project is archived at:

https://doi.org/10.17632/r97b37bfjg.1

## Environment

Python 3.12 is recommended. The archived package records the following exact core dependencies:

```text
numpy==2.3.5
scipy==1.17.0
pytest==9.1.1
```

## Canonical workflow

1. Download and extract the Version 1 archive from Mendeley Data.
2. Create an isolated Python 3.12 environment.
3. Install the dependencies from the archive's `requirements.txt`.
4. Run the package-integrity verification script.
5. Run the complete reproduction wrapper.
6. Review the generated run summary and task logs.

The archived workflow is designed to run without network access after dependencies are installed. Reproduction begins from included provenance-labelled transcribed rows; copyrighted source PDFs are not redistributed.

## Interpretation boundary

The archive includes research-case records, numerical regression checks, method audits, and source assessments. These categories are not interchangeable.

In particular:

- numerical regression checks are not independent scientific validations;
- preserved historical tests establish software behavior, not experimental confirmation;
- synthetic method audits use known generators and are not blind hold-out experiments;
- no new experimental measurements are introduced by the archive.

The Mendeley Data version is the canonical release if this GitHub companion and the archive ever differ.
