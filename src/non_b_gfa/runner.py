"""Top-level orchestration for parsing input, running detectors, and writing outputs."""

from __future__ import annotations

import subprocess
from pathlib import Path

from .defaults import RunConfig
from .detectors import SequenceContext, complement, find_apr, find_dr, find_gq, find_ir, find_mr, find_str, find_zdna, is_subset, reverse_complement
from .fasta import read_fasta_records
from .models import Repeat
from .writers import write_gff, write_tsv

SUFFIXES = {
    "I": "_IR",
    "M": "_MR",
    "D": "_DR",
    "G": "_GQ",
    "Z": "_Z",
    "S": "_STR",
    "A": "_APR",
}


def _seq_title(record_title: str, chrom: str | None, fasta_index: int, fasta_count: int) -> str:
    """Derive the output sequence name for one FASTA record."""

    if chrom:
        title = chrom
        if fasta_count > 1:
            title = f"{title}_{fasta_index}"
        return title
    return record_title.split()[0]


def _open_outputs(config: RunConfig) -> dict[str, tuple[Path, Path, object, object]]:
    outputs: dict[str, tuple[Path, Path, object, object]] = {}
    for motif in config.enabled_motifs():
        gff_path = Path(f"{config.output_prefix}{SUFFIXES[motif]}.gff")
        tsv_path = Path(f"{config.output_prefix}{SUFFIXES[motif]}.tsv")
        outputs[motif] = (
            gff_path,
            tsv_path,
            gff_path.open("w", encoding="ascii"),
            tsv_path.open("w", encoding="ascii"),
        )
    return outputs


def _run_detectors(config: RunConfig, context: SequenceContext) -> dict[str, list[Repeat]]:
    """Run all enabled motif detectors and subset annotations."""

    results: dict[str, list[Repeat]] = {}
    if config.findIR:
        repeats = find_ir(context, config.minIRrep, config.maxIRspacer, config.shortIRcut, config.shortIRspacer)
        if config.findCruciform:
            is_subset(repeats, "I", config.maxCruciformSpacer, config.minCruciformRep, context.dna)
        results["I"] = repeats
    if config.findGQ:
        results["G"] = find_gq(context, config.minGQrep, config.maxGQspacer)
    if config.findMR:
        repeats = find_mr(context, config.minMRrep, config.maxMRspacer)
        if config.findTriplex:
            is_subset(repeats, "M", config.maxTriplexSpacer, config.minTriplexYRpercent, context.dna)
        results["M"] = repeats
    if config.findDR:
        repeats = find_dr(context, config.minDRrep, config.maxDRrep, config.maxDRspacer)
        if config.findSlipped:
            is_subset(repeats, "D", config.maxSlippedSpacer, -999, context.dna)
        results["D"] = repeats
    if config.findZ:
        repeats = find_zdna(context, config.minZlen)
        if config.findKVzdna:
            is_subset(repeats, "Z", -999, config.minKVscore, context.dna)
        results["Z"] = repeats
    if config.findSTR:
        results["S"] = find_str(context, config.minSTR, config.maxSTR, config.minSTRbp, config.minSTRreps)
    if config.findAPR:
        results["A"] = find_apr(context, config.minAPRlen, config.maxAPRlen, config.minATracts)
    return results


def run_analysis(config: RunConfig) -> None:
    """Execute one configured analysis run and write all enabled outputs."""

    records = read_fasta_records(config.sequence_file)
    outputs = _open_outputs(config)
    counters = {motif: 0 for motif in outputs}
    try:
        for fasta_index, record in enumerate(records, start=1):
            context = SequenceContext(
                dna=record.sequence,
                dna2=reverse_complement(record.sequence),
                dna3=complement(record.sequence),
            )
            seq_title = _seq_title(record.title, config.chrom, fasta_index, len(records))
            results = _run_detectors(config, context)
            for motif, repeats in results.items():
                if not repeats:
                    continue
                _, _, gff_handle, tsv_handle = outputs[motif]
                write_gff(gff_handle, repeats, seq_title, motif, context.dna, context.dna2)
                counters[motif] += 1
                write_tsv(tsv_handle, repeats, seq_title, motif, context.dna, context.dna2, write_header=counters[motif] == 1)
    finally:
        for _, _, gff_handle, tsv_handle in outputs.values():
            gff_handle.close()
            tsv_handle.close()

    if config.do_wget:
        subprocess.run(
            [
                "wget",
                "-S",
                "-",
                f"http://nonb.abcc.ncifcrf.gov/modules/nBMSTc/controllers/notify.php?file={config.output_prefix}",
            ],
            check=False,
        )
    if config.do_chmod:
        for motif, (gff_path, tsv_path, _, _) in outputs.items():
            del motif
            gff_path.chmod(0o664)
            tsv_path.chmod(0o664)
