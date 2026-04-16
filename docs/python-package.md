# Python package guide

## Overview

`non_b_gfa` is the Python reimplementation of the legacy `gfa` motif finder. It keeps the original command-line behavior while exposing a small Python API for programmatic use.

## Public API

### `RunConfig`

`RunConfig` defines the input file, output prefix, search thresholds, skip flags, and optional side effects such as `wget` notification and `chmod`.

Minimal example:

```python
from pathlib import Path

from non_b_gfa import RunConfig

config = RunConfig(
    sequence_file=Path("input.fasta"),
    output_prefix=Path("results/run1"),
    do_wget=False,
)
```

### `run_analysis`

`run_analysis(config)` reads one or more FASTA records, runs the enabled motif detectors, and writes motif-specific `.gff` and `.tsv` output files using the configured prefix.

```python
from pathlib import Path

from non_b_gfa import RunConfig, run_analysis

run_analysis(
    RunConfig(
        sequence_file=Path("input.fasta"),
        output_prefix=Path("output_prefix"),
        do_wget=False,
    )
)
```

## Command-line usage

After installation, the package exposes a `gfa` command:

```bash
gfa -seq input.fasta -out results/output_prefix
```

The CLI accepts the same switch names as the legacy C program, including motif thresholds and `-skip*` flags. The legacy completion callback is disabled by default in the Python implementation; pass `-doWGET` if you need to opt in.

## Output behavior

- each enabled motif writes one `.gff` file and one `.tsv` file
- output suffixes match the historical tool, such as `_IR`, `_MR`, `_DR`, `_GQ`, `_Z`, `_STR`, and `_APR`
- when multiple FASTA records are present, output is appended into the same motif files and record names are derived from the FASTA title or the supplied `-chrom` value

## Internal modules

- `non_b_gfa.fasta` parses FASTA files into `FastaRecord` objects
- `non_b_gfa.detectors` implements the motif-finding algorithms and subset tagging
- `non_b_gfa.writers` renders GFF and TSV records
- `non_b_gfa.runner` coordinates parsing, detection, and writing

## Validation

The Python implementation is regression-tested against the bundled `test_files.tar` sample outputs.
