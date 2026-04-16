"""FASTA parsing helpers."""

from __future__ import annotations

from pathlib import Path

from .models import FastaRecord


def read_fasta_records(path: Path) -> list[FastaRecord]:
    """Read a FASTA file into normalized lowercase records.

    Non-alphabetic characters inside sequence lines are ignored to match the
    permissive behavior of the legacy tool.
    """

    records: list[FastaRecord] = []
    title: str | None = None
    sequence_parts: list[str] = []

    with path.open("r", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if title is not None:
                    records.append(FastaRecord(title=title, sequence="".join(sequence_parts).lower()))
                title = line[1:].strip()
                sequence_parts = []
                continue
            sequence_parts.append("".join(ch for ch in line if ch.isalpha()))

    if title is not None:
        records.append(FastaRecord(title=title, sequence="".join(sequence_parts).lower()))

    if not records or path.read_text(encoding="ascii").lstrip()[:1] != ">":
        raise ValueError("Sequence file is empty or not in FASTA format.")

    return records
