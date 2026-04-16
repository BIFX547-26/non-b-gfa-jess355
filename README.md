# non-b-gfa

## Overview

`non-b-gfa` detects non-B DNA motifs from FASTA input and writes motif annotations as GFF and TSV output files. The project preserves the behavior and output format of the original motif-finding tool while providing a cleaner Python implementation. It is a Python refactor of a legacy codebase that was originally written in C.

## Project Structure

- `src/non_b_gfa/` - main Python package containing the CLI, FASTA parsing, motif detection logic, configuration, and output writers
- `tests/` - unit and regression tests for the Python implementation
- `data/` - sample FASTA input and reference output files used for comparison and validation
- `legacy_c/` - original C source code kept for reference during the refactor

## Installation

Install the project in editable mode with:

```bash
pip install -e .
```

## Usage

Run the CLI with:

```bash
PYTHONPATH=src python3 -m non_b_gfa.cli -seq data/gfa_test.fasta -out test -chrom test
```

This will generate multiple `.gff` and `.tsv` files prefixed with the value passed to `-out`.

## Output

The program generates motif-specific `.gff` and `.tsv` files based on the input sequence. Different motif types produce separate output files, such as inverted repeats, mirror repeats, direct repeats, G-quadruplexes, Z-DNA, short tandem repeats, and A-phased repeats.

## Validation

The Python implementation was validated by running the program on the provided `gfa_test.fasta` input and comparing all generated `.gff` and `.tsv` files against the reference outputs in the `data/` directory using `diff`.

The outputs matched the reference results, with only minor capitalization differences. This confirms that the refactored Python version preserves the behavior of the original implementation.

## Notes

The original implementation included functionality to submit jobs to an external web service using `wget` and `notify.php`. This behavior was removed in the Python refactor.

All processing in the Python version is performed locally and does not rely on external services.
