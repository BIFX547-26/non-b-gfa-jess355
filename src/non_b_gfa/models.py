"""Shared dataclasses used across parsing, detection, and writing."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class Repeat:
    """A detected motif record.

    The fields mirror the layout used by the legacy C implementation so the
    Python port can preserve its output semantics.
    """

    start: int = 0
    loop: int = 0
    len: int = 0
    num: int = 0
    end: int = 0
    sub: int = 0
    strand: int = 0
    special: int = 0


@dataclass(slots=True)
class GIsland:
    """A consecutive run of G or C bases used during G-quadruplex detection."""

    strt: int = 0
    len: int = 0


@dataclass(slots=True)
class PotentialBentDNA:
    """Intermediate A-tract information used while building APR calls."""

    a_center: float = 0.0
    strt: int = 0
    end: int = 0


@dataclass(slots=True)
class FastaRecord:
    """A single FASTA title/sequence pair."""

    title: str
    sequence: str
