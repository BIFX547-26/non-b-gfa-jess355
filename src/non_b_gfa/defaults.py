"""Runtime configuration defaults for the Python package."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class RunConfig:
    """Configuration for a single analysis run."""

    sequence_file: Path
    output_prefix: Path
    chrom: str | None = None
    minGQrep: int = 3
    maxGQspacer: int = 7
    minMRrep: int = 10
    maxMRspacer: int = 100
    minIRrep: int = 6
    maxIRspacer: int = 100
    shortIRcut: int = 9
    shortIRspacer: int = 4
    minDRrep: int = 10
    maxDRrep: int = 300
    maxDRspacer: int = 10
    minATracts: int = 3
    minATractSep: int = 10
    maxATractSep: int = 11
    maxAPRlen: int = 9
    minAPRlen: int = 3
    minZlen: int = 10
    minSTR: int = 1
    maxSTR: int = 9
    minSTRbp: int = 10
    minSTRreps: int = 3
    minCruciformRep: int = 10
    maxCruciformSpacer: int = 4
    minTriplexYRpercent: int = 10
    maxTriplexSpacer: int = 8
    maxSlippedSpacer: int = 0
    minKVscore: int = 33
    findAPR: bool = True
    findSTR: bool = True
    findDR: bool = True
    findMR: bool = True
    findIR: bool = True
    findGQ: bool = True
    findZ: bool = True
    findSlipped: bool = True
    findCruciform: bool = True
    findTriplex: bool = True
    findKVzdna: bool = True
    do_wget: bool = False
    do_chmod: bool = False

    def enabled_motifs(self) -> list[str]:
        """Return motif codes enabled for the current run."""

        motifs: list[str] = []
        if self.findIR:
            motifs.append("I")
        if self.findMR:
            motifs.append("M")
        if self.findDR:
            motifs.append("D")
        if self.findGQ:
            motifs.append("G")
        if self.findZ:
            motifs.append("Z")
        if self.findSTR:
            motifs.append("S")
        if self.findAPR:
            motifs.append("A")
        return motifs
