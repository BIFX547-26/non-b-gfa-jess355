"""Writers for the legacy GFF and TSV output formats."""

from __future__ import annotations

from typing import TextIO

from .models import Repeat


def _composition(segment: str) -> tuple[int, int, int, int]:
    return segment.count("a"), segment.count("c"), segment.count("g"), segment.count("t")


def _get_seq_fields(repeat: Repeat, motif: str, dna: str, dna2: str, total_bases: int) -> tuple[str, tuple[int, int, int, int]]:
    if motif in {"D", "M", "I", "S"}:
        start = repeat.start - 1
        stop = start + repeat.len
    else:
        start = repeat.start - 1
        stop = repeat.end - 1

    if repeat.strand == 0:
        sequence = dna[repeat.start - 1 : repeat.end]
        comp_segment = dna[start:stop]
    else:
        start = total_bases - start - 1
        stop = total_bases - stop - 1
        sequence = dna2[stop : stop + (repeat.end - repeat.start + 1)]
        comp_segment = dna2[stop : start + 1]
    return sequence, _composition(comp_segment)


def _rep_names(motif: str) -> tuple[str, str]:
    names = {
        "M": ("Mirror_Repeat", "MR"),
        "I": ("Inverted_Repeat", "IR"),
        "D": ("Direct_Repeat", "DR"),
        "G": ("G_Quadruplex_Motif", "GQ"),
        "Z": ("Z_DNA_Motif", "ZDNA"),
        "A": ("A_Phased_Repeat", "APR"),
        "S": ("Short_Tandem_Repeat", "STR"),
    }
    return names[motif]


def write_gff(handle: TextIO, repeats: list[Repeat], chrom: str, motif: str, dna: str, dna2: str) -> None:
    """Write detected repeats in the historical GFF output format."""

    rep_type, short_type = _rep_names(motif)
    total_bases = len(dna)
    handle.write("##gff-version 3\n")
    for repeat in repeats:
        sequence, (a_count, c_count, g_count, t_count) = _get_seq_fields(repeat, motif, dna, dna2, total_bases)
        fields = [chrom, "ABCC", rep_type, str(repeat.start), str(repeat.end), "."]
        if motif == "G":
            fields.extend(["+" if repeat.strand == 0 else "-", "."])
            attrs = [
                f"ID={chrom}_{repeat.start}_{repeat.end}_{short_type}",
                f"islands={repeat.sub}",
                f"runs={repeat.num}",
                f"max={repeat.len}",
                f"composition={a_count}A/{c_count}C/{g_count}G/{t_count}T",
                f"sequence={sequence}",
            ]
        else:
            fields.extend(["+", "."])
            attrs = [f"ID={chrom}_{repeat.start}_{repeat.end}_{short_type}"]
            if motif == "Z":
                attrs.extend([f"length={repeat.len}", f"score={repeat.loop}"])
            elif motif == "A":
                attrs.append(f"tracts={repeat.num}")
            elif motif == "S":
                attrs.extend([f"length={repeat.len}", f"x{repeat.num}+{repeat.sub}", f"type={repeat.loop}"])
            elif motif == "D":
                attrs.extend([f"spacer={repeat.loop}", f"repeat={repeat.len}", f"x{repeat.num}+{repeat.sub}"])
            else:
                attrs.extend([f"spacer={repeat.loop}", f"repeat={repeat.len}", f"perms={repeat.num}", f"minloop={repeat.sub}"])
            attrs.extend([f"composition={a_count}A/{c_count}C/{g_count}G/{t_count}T", f"sequence={sequence}"])
            if motif in {"D", "M", "I", "Z"}:
                attrs.append(f"subset={repeat.special}")
        handle.write("\t".join(fields + [";".join(attrs)]) + "\n")


def write_tsv(handle: TextIO, repeats: list[Repeat], chrom: str, motif: str, dna: str, dna2: str, write_header: bool) -> None:
    """Write detected repeats in the historical TSV output format."""

    rep_type, _ = _rep_names(motif)
    total_bases = len(dna)
    if write_header:
        handle.write("Sequence_name\tSource\tType\tStart\tStop\tLength\tScore\tStrand\tRepeat\tSpacer\t")
        detail_header = {
            "D": "Repeated",
            "M": "Permutations",
            "I": "Permutations",
            "G": "nIslands/nRuns/maxGQ",
            "Z": "KVScore",
            "A": "Tracts",
            "S": "Repeated",
        }[motif]
        handle.write(f"{detail_header}\tSubset\tComposition\tSequence\n")

    for repeat in repeats:
        sequence, (a_count, c_count, g_count, t_count) = _get_seq_fields(repeat, motif, dna, dna2, total_bases)
        spacer = "NA" if motif in {"G", "Z"} else str(repeat.loop)
        detail = {
            "D": f"X{repeat.num}+{repeat.sub}",
            "M": str(repeat.num),
            "G": f"{repeat.sub}I/{repeat.num}R/{repeat.len}M",
            "I": str(repeat.num),
            "Z": str(repeat.loop),
            "A": str(repeat.num),
            "S": f"X{repeat.num}+{repeat.sub}",
        }[motif]
        subset = str(repeat.special) if motif in {"D", "M", "I", "Z"} else "NA"
        handle.write(
            "\t".join(
                [
                    chrom,
                    "ABCC",
                    rep_type,
                    str(repeat.start),
                    str(repeat.end),
                    str(repeat.end - repeat.start + 1),
                    "NA",
                    "+" if repeat.strand == 0 else "-",
                    str(repeat.len),
                    spacer,
                    detail,
                    subset,
                    f"{a_count}A/{c_count}C/{g_count}G/{t_count}T",
                    sequence,
                ]
            )
            + "\n"
        )
