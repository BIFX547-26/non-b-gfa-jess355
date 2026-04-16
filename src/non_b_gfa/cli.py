"""Command-line interface for the Python `gfa` package."""

from __future__ import annotations

import argparse
from pathlib import Path

from .defaults import RunConfig
from .runner import run_analysis


def build_parser() -> argparse.ArgumentParser:
    """Build the legacy-compatible CLI parser."""

    parser = argparse.ArgumentParser(prog="gfa")
    parser.add_argument("-seq", required=True, type=Path, dest="sequence_file")
    parser.add_argument("-out", required=True, type=Path, dest="output_prefix")
    parser.add_argument("-chrom")
    parser.add_argument("-minGQrep", type=int, default=3)
    parser.add_argument("-maxGQspacer", type=int, default=7)
    parser.add_argument("-minMRrep", type=int, default=10)
    parser.add_argument("-maxMRspacer", type=int, default=100)
    parser.add_argument("-minIRrep", type=int, default=6)
    parser.add_argument("-maxIRspacer", type=int, default=100)
    parser.add_argument("-shortIRcut", type=int, default=9)
    parser.add_argument("-shortIRspacer", type=int, default=4)
    parser.add_argument("-minDRrep", type=int, default=10)
    parser.add_argument("-maxDRrep", type=int, default=300)
    parser.add_argument("-maxDRspacer", type=int, default=10)
    parser.add_argument("-minATracts", type=int, default=3)
    parser.add_argument("-minATractSep", type=int, default=10)
    parser.add_argument("-maxATractSep", type=int, default=11)
    parser.add_argument("-maxAPRlen", type=int, default=9)
    parser.add_argument("-minAPRlen", type=int, default=3)
    parser.add_argument("-minZlen", type=int, default=10)
    parser.add_argument("-minSTR", type=int, default=1)
    parser.add_argument("-maxSTR", type=int, default=9)
    parser.add_argument("-minSTRbp", type=int, default=10)
    parser.add_argument("-minCruciformRep", type=int, default=10)
    parser.add_argument("-maxCruciformSpacer", type=int, default=4)
    parser.add_argument("-minTriplexYRpercent", type=int, default=10)
    parser.add_argument("-maxTriplexSpacer", type=int, default=8)
    parser.add_argument("-maxSlippedSpacer", type=int, default=0)
    parser.add_argument("-skipAPR", action="store_true")
    parser.add_argument("-skipSTR", action="store_true")
    parser.add_argument("-skipDR", action="store_true")
    parser.add_argument("-skipMR", action="store_true")
    parser.add_argument("-skipIR", action="store_true")
    parser.add_argument("-skipGQ", action="store_true")
    parser.add_argument("-skipZ", action="store_true")
    parser.add_argument("-skipSlipped", action="store_true")
    parser.add_argument("-skipCruciform", action="store_true")
    parser.add_argument("-skipTriplex", action="store_true")
    parser.add_argument("-skipKVzdna", action="store_true")
    parser.add_argument("-skipWGET", action="store_true")
    parser.add_argument("-doWGET", action="store_true")
    parser.add_argument("-doCHMOD", action="store_true")
    return parser


def main() -> None:
    """Parse CLI arguments and execute the analysis run."""

    args = build_parser().parse_args()
    config = RunConfig(
        sequence_file=args.sequence_file,
        output_prefix=args.output_prefix,
        chrom=args.chrom,
        minGQrep=args.minGQrep,
        maxGQspacer=args.maxGQspacer,
        minMRrep=args.minMRrep,
        maxMRspacer=args.maxMRspacer,
        minIRrep=args.minIRrep,
        maxIRspacer=args.maxIRspacer,
        shortIRcut=args.shortIRcut,
        shortIRspacer=args.shortIRspacer,
        minDRrep=args.minDRrep,
        maxDRrep=args.maxDRrep,
        maxDRspacer=args.maxDRspacer,
        minATracts=args.minATracts,
        minATractSep=args.minATractSep,
        maxATractSep=args.maxATractSep,
        maxAPRlen=args.maxAPRlen,
        minAPRlen=args.minAPRlen,
        minZlen=args.minZlen,
        minSTR=args.minSTR,
        maxSTR=args.maxSTR,
        minSTRbp=args.minSTRbp,
        minCruciformRep=args.minCruciformRep,
        maxCruciformSpacer=args.maxCruciformSpacer,
        minTriplexYRpercent=args.minTriplexYRpercent,
        maxTriplexSpacer=args.maxTriplexSpacer,
        maxSlippedSpacer=args.maxSlippedSpacer,
        findAPR=not args.skipAPR,
        findSTR=not args.skipSTR,
        findDR=not args.skipDR,
        findMR=not args.skipMR,
        findIR=not args.skipIR,
        findGQ=not args.skipGQ,
        findZ=not args.skipZ,
        findSlipped=not args.skipSlipped,
        findCruciform=not args.skipCruciform,
        findTriplex=not args.skipTriplex,
        findKVzdna=not args.skipKVzdna,
        do_wget=args.doWGET and not args.skipWGET,
        do_chmod=args.doCHMOD,
    )
    run_analysis(config)


if __name__ == "__main__":
    main()
