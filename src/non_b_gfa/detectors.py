"""Motif detection and subset-tagging logic ported from the legacy C implementation."""

from __future__ import annotations

from dataclasses import dataclass
from math import floor

from .models import GIsland, PotentialBentDNA, Repeat


def complement(sequence: str) -> str:
    """Return the base complement of a lowercase DNA sequence."""

    table = str.maketrans("acgtn", "tgcan")
    return sequence.translate(table)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a lowercase DNA sequence."""

    return complement(sequence)[::-1]


def _base(sequence: str, index: int, default: str = "\0") -> str:
    if 0 <= index < len(sequence):
        return sequence[index]
    return default


@dataclass(slots=True)
class SequenceContext:
    """Sequence buffers shared by the motif detectors."""

    dna: str
    dna2: str
    dna3: str

    @property
    def total_bases(self) -> int:
        return len(self.dna)


def _shift_delete(repeats: list[Repeat], to_remove: int) -> None:
    del repeats[to_remove]


def find_ir(context: SequenceContext, mincrf: int, cspacer: int, cut: int, short_spacer: int) -> list[Repeat]:
    """Find inverted repeats."""

    dna = context.dna
    dna3 = context.dna3
    total_bases = context.total_bases
    repeats: list[Repeat] = []

    for strti in range(mincrf, total_bases - mincrf + 1):
        max_sp = min(cspacer, total_bases - (strti + mincrf))
        right_shifted = False
        left_shifted = False
        for sp in range(max_sp + 1):
            i = strti
            j = strti + sp + 1
            k = 0
            while j < total_bases and i >= 0 and _base(dna3, j) == _base(dna, i) and _base(dna, j) != "n":
                k += 1
                j += 1
                i -= 1
            if k < mincrf:
                continue
            if k <= cut and sp > short_spacer:
                continue

            tmp_start = (strti - k) + 2
            tmp_stop = strti + k + sp + 1
            if not repeats:
                repeats.append(Repeat(start=tmp_start, sub=tmp_stop, len=k, loop=sp, num=1, end=tmp_stop, strand=0))
                continue

            while repeats and repeats[-1].end <= tmp_stop and repeats[-1].start >= tmp_start and repeats[-1].len < k:
                repeats.pop()
                right_shifted = False
                left_shifted = False
            while repeats and repeats[-1].end >= tmp_stop and repeats[-1].start <= tmp_start and repeats[-1].len < k:
                repeats.pop()
                right_shifted = False
                left_shifted = False
            if repeats and repeats[-1].end <= tmp_stop and repeats[-1].start >= tmp_start and repeats[-1].len > k:
                right_shifted = False
                left_shifted = False
                continue
            if repeats and repeats[-1].end >= tmp_stop and repeats[-1].start <= tmp_start and repeats[-1].len > k:
                right_shifted = False
                left_shifted = False
                continue
            if repeats and tmp_stop == repeats[-1].end and k == repeats[-1].len and not right_shifted:
                left_shifted = True
                right_shifted = False
                repeats[-1].num += 1
                repeats[-1].sub = tmp_start
                continue
            if repeats and tmp_start == repeats[-1].start and k == repeats[-1].len and not left_shifted:
                right_shifted = True
                left_shifted = False
                if len(repeats) >= 2 and repeats[-2].end <= tmp_stop and repeats[-2].start >= tmp_start and repeats[-2].len < k:
                    repeats[-2] = repeats[-1]
                    repeats.pop()
                repeats[-1].num += 1
                repeats[-1].end = tmp_stop
                continue

            right_shifted = False
            left_shifted = False
            repeats.append(Repeat(start=tmp_start, sub=tmp_stop, len=k, loop=sp, num=1, end=tmp_stop, strand=0))
            c_back = 1
            while c_back <= 10 and c_back < len(repeats):
                prev = repeats[-1 - c_back]
                curr = repeats[-1]
                included = ((prev.end >= curr.end and prev.start <= curr.start) or (prev.end <= curr.end and prev.start >= curr.start))
                if not included:
                    c_back += 1
                    continue
                if prev.len == curr.len:
                    if prev.loop > curr.loop:
                        _shift_delete(repeats, len(repeats) - 1 - c_back)
                        curr = repeats[-1]
                elif prev.len < curr.len:
                    _shift_delete(repeats, len(repeats) - 1 - c_back)
                    curr = repeats[-1]
                else:
                    c_back += 1
                    continue
                repeats.pop()
                right_shifted = False
                left_shifted = False
                c_back = 1
    return repeats


def find_mr(context: SequenceContext, minmir: int, mspacer: int) -> list[Repeat]:
    """Find mirror repeats."""

    dna = context.dna
    total_bases = context.total_bases
    repeats: list[Repeat] = []

    for strti in range(minmir, total_bases - minmir + 1):
        max_sp = min(mspacer, total_bases - (strti + minmir))
        right_shifted = False
        left_shifted = False
        for sp in range(max_sp + 1):
            i = strti
            j = strti + sp + 1
            k = 0
            while j < total_bases and i >= 0 and _base(dna, i) == _base(dna, j) and _base(dna, j) != "n":
                k += 1
                j += 1
                i -= 1
            if k < minmir:
                continue

            tmp_start = (strti - k) + 2
            tmp_stop = strti + k + sp + 1
            if not repeats:
                repeats.append(Repeat(start=tmp_start, sub=tmp_stop, len=k, loop=sp, num=1, end=tmp_stop, strand=0))
                continue

            while repeats and repeats[-1].end <= tmp_stop and repeats[-1].start >= tmp_start and repeats[-1].len < k:
                repeats.pop()
                right_shifted = False
                left_shifted = False
            while repeats and repeats[-1].end >= tmp_stop and repeats[-1].start <= tmp_start and repeats[-1].len < k:
                repeats.pop()
                right_shifted = False
                left_shifted = False
            if repeats and repeats[-1].end <= tmp_stop and repeats[-1].start >= tmp_start and repeats[-1].len > k:
                right_shifted = False
                left_shifted = False
                continue
            if repeats and repeats[-1].end >= tmp_stop and repeats[-1].start <= tmp_start and repeats[-1].len > k:
                right_shifted = False
                left_shifted = False
                continue
            if repeats and tmp_stop == repeats[-1].end and k == repeats[-1].len and not right_shifted:
                left_shifted = True
                right_shifted = False
                repeats[-1].num += 1
                repeats[-1].sub = tmp_start
                continue
            if repeats and tmp_start == repeats[-1].start and k == repeats[-1].len and not left_shifted:
                right_shifted = True
                left_shifted = False
                if len(repeats) >= 2 and repeats[-2].end <= tmp_stop and repeats[-2].start >= tmp_start and repeats[-2].len < k:
                    repeats[-2] = repeats[-1]
                    repeats.pop()
                repeats[-1].num += 1
                repeats[-1].end = tmp_stop
                continue

            right_shifted = False
            left_shifted = False
            repeats.append(Repeat(start=tmp_start, sub=tmp_stop, len=k, loop=sp, num=1, end=tmp_stop, strand=0))
            c_back = 1
            while c_back <= 5 and c_back < len(repeats):
                prev = repeats[-1 - c_back]
                curr = repeats[-1]
                included = ((prev.end >= curr.end and prev.start <= curr.start) or (prev.end <= curr.end and prev.start >= curr.start))
                if not included:
                    c_back += 1
                    continue
                if prev.len == curr.len:
                    if prev.loop > curr.loop:
                        _shift_delete(repeats, len(repeats) - 1 - c_back)
                        curr = repeats[-1]
                elif prev.len < curr.len:
                    _shift_delete(repeats, len(repeats) - 1 - c_back)
                    curr = repeats[-1]
                else:
                    c_back += 1
                    continue
                repeats.pop()
                right_shifted = False
                left_shifted = False
                c_back = 1
    return repeats


def find_dr(context: SequenceContext, mindir: int, maxdir: int, dspacer: int) -> list[Repeat]:
    """Find direct repeats."""

    dna = context.dna
    total_bases = context.total_bases
    repeats: list[Repeat] = []
    end = 0
    lasti = total_bases - (mindir * 2)

    strti = 0
    while strti <= lasti:
        while strti <= lasti and _base(dna, strti) == "n":
            strti += 1
        if strti >= lasti:
            break
        for size in range(maxdir, mindir - 1, -1):
            if ((size * 2) + dspacer) <= (end - strti):
                continue
            sp_min = max(0, ((end - strti) - (size * 2)) + 2)
            sp_max = min(dspacer, lasti - strti)
            for sp in range(sp_min, sp_max + 1):
                j = strti + size + sp
                i = strti
                k = 0
                while j < total_bases and _base(dna, i) == _base(dna, j) and k < size and _base(dna, i) != "n":
                    k += 1
                    j += 1
                    i += 1
                if k != size:
                    continue
                totlen = k
                if sp == 0:
                    while j < total_bases and i < total_bases and _base(dna, i) == _base(dna, j):
                        totlen += 1
                        j += 1
                        i += 1
                repeats.append(
                    Repeat(
                        start=strti + 1,
                        len=size,
                        loop=sp,
                        num=totlen // size,
                        end=j,
                        sub=totlen % size,
                        strand=0,
                    )
                )
                end = j - 1
                break
            else:
                continue
            break
        strti += 1
    return repeats


def get_gislands(context: SequenceContext, min_gq: int) -> tuple[list[GIsland], list[GIsland]]:
    """Collect G-islands on the forward strand and C-islands on the reverse strand."""

    dna = context.dna
    gisle: list[GIsland] = []
    rcgisle: list[GIsland] = []
    ngs = 0
    ncs = 0

    for i in range(context.total_bases + 1):
        base = _base(dna, i)
        if base == "g":
            ngs += 1
        else:
            if ngs >= min_gq:
                gisle.append(GIsland(strt=i - ngs + 1, len=ngs))
            ngs = 0
        if base == "c":
            ncs += 1
        else:
            if ncs >= min_gq:
                rcgisle.append(GIsland(strt=i - ncs + 1, len=ncs))
            ncs = 0
    return gisle, rcgisle


def find_gq(context: SequenceContext, min_gq: int, max_gq_spacer: int) -> list[Repeat]:
    """Find G-quadruplex motifs."""

    repeats: list[Repeat] = []
    for strand, islands in enumerate(get_gislands(context, min_gq)):
        n_ils = len(islands)
        i = 0
        while i < n_ils:
            con_ils = 1
            npos = int(floor((islands[i].len + 1) / (min_gq + 1)))
            i2 = i + 1
            while i2 < n_ils and (islands[i2].strt - (islands[i2 - 1].strt + islands[i2 - 1].len)) <= max_gq_spacer:
                con_ils += 1
                npos += int(floor((islands[i2].len + 1) / (min_gq + 1)))
                i2 += 1
            if npos >= 4:
                max_gq = min_gq
                for j in range(i, i2):
                    for k in range(islands[j].len, max_gq, -1):
                        npos_max = int(floor((islands[j].len + 1) / (k + 1)))
                        for m in range(j + 1, i2):
                            npos_max += int(floor((islands[m].len + 1) / (k + 1)))
                            if npos_max >= 4:
                                max_gq = k
                                break
                            if int(floor((islands[m].len + 1) / (k + 1))) == 0:
                                if m + 1 < i2 and islands[m + 1].strt > (islands[m - 1].strt + islands[m - 1].len + max_gq_spacer):
                                    break
                        else:
                            continue
                        break
                repeats.append(
                    Repeat(
                        start=islands[i].strt,
                        num=npos,
                        sub=con_ils,
                        len=max_gq,
                        end=(islands[i2 - 1].strt + islands[i2 - 1].len) - 1,
                        strand=strand,
                    )
                )
            i += con_ils
    return repeats


def _pupy(sequence: str, pos: int) -> int:
    first = _base(sequence, pos)
    second = _base(sequence, pos + 1)
    if first == "a" and second == "c":
        return 3
    if first == "t" and second == "g":
        return 3
    if first == "c":
        if second == "g":
            return 25
        if second == "a":
            return 3
    if first == "g":
        if second == "c":
            return 25
        if second == "t":
            return 3
    return 0


def find_zdna(context: SequenceContext, min_z: int) -> list[Repeat]:
    """Find Z-DNA candidate runs."""

    dna = context.dna
    total_bases = context.total_bases
    repeats: list[Repeat] = []
    i = 0
    npy = 1
    kvsum = 0

    while i < (total_bases - min_z):
        tmp_ppy = _pupy(dna, i)
        if tmp_ppy > 0:
            npy += 1
            kvsum += tmp_ppy
        else:
            if npy >= min_z:
                repeats.append(Repeat(start=i - npy + 2, len=npy, loop=kvsum // 2, num=0, end=i + 1, sub=0, strand=0))
            npy = 1
            kvsum = 0
        i += 1
    return repeats


def _non_b_str(sequence: str, start: int, length: int) -> int:
    code = 0
    is_even = length % 2 == 0
    is_symmetric = True
    is_pupy = True
    is_comp = True
    j = start + length - 2

    if length >= 2:
        for i in range(0, (length // 2)):
            left = _base(sequence, start + i - 1)
            right = _base(sequence, j)
            if left != right:
                is_symmetric = False
            if is_even:
                comp_map = {"a": "t", "t": "a", "c": "g", "g": "c"}
                if comp_map.get(left) != right:
                    is_comp = False
            else:
                is_comp = False
            j -= 1
        for i in range(start, length + start - 1):
            curr = _base(sequence, i)
            prev = _base(sequence, i - 1)
            if curr in {"a", "g"} and prev in {"a", "g"}:
                is_pupy = False
            if curr in {"t", "c"} and prev in {"t", "c"}:
                is_pupy = False
    else:
        is_comp = False
        is_symmetric = True
        is_pupy = False

    if is_even:
        code += 1
    if is_pupy:
        code += 2
    if is_symmetric:
        code += 4
    if is_comp:
        code += 8
    return code


def find_str(context: SequenceContext, min_str: int, max_str: int, min_strlen: int, min_reps: int) -> list[Repeat]:
    """Find short tandem repeats."""

    dna = context.dna
    total_bases = context.total_bases
    repeats: list[Repeat] = []
    i = 0
    while i < (total_bases - min_strlen):
        while i != (total_bases - 1) and _base(dna, i) == "n":
            i += 1
        found = False
        for rpsz in range(min_str, max_str + 1):
            reps = 1
            j = i + rpsz
            while j + rpsz < total_bases and dna[i : i + rpsz] == dna[j : j + rpsz]:
                reps += 1
                j += rpsz
            if reps < min_reps:
                continue
            remainder = 0
            rs = i
            re = j
            while rs < total_bases and re < total_bases and _base(dna, rs) == _base(dna, re):
                remainder += 1
                rs += 1
                re += 1
            if ((reps * rpsz) + remainder) < min_strlen:
                continue
            if repeats and repeats[-1].end >= re:
                continue
            repeats.append(
                Repeat(
                    start=i + 1,
                    end=re,
                    num=reps,
                    loop=_non_b_str(dna, i + 1, rpsz),
                    len=rpsz,
                    sub=remainder,
                    strand=0,
                )
            )
            i = re - min_strlen + 1
            found = True
            break
        i += 1 if not found else 0
    return repeats


def _get_atracts(context: SequenceContext, min_at: int, max_at: int) -> list[PotentialBentDNA]:
    dna = context.dna
    dna2 = context.dna2
    total_bases = context.total_bases
    patrs: list[PotentialBentDNA] = []
    i = 0
    n_as = 0

    while i < total_bases:
        if _base(dna, i) in {"a", "t"}:
            n_as += 1
        else:
            if min_at <= n_as <= max_at:
                strt = i - n_as + 1
                at_end = strt + n_as

                max_atlen = max_tlen = alen = tlen = atlen = talen = max_atend = 0
                max_atlen_rc = max_tlen_rc = alen_rc = tlen_rc = atlen_rc = talen_rc = max_atend_rc = 0

                n_rc = total_bases - at_end
                for n in range(strt - 1, at_end - 1):
                    n_rc += 1
                    if _base(dna, n) == "a":
                        tlen = 0
                        talen = 0
                        if _base(dna, n - 1) == "t":
                            alen = 0
                            atlen = 0
                        else:
                            alen += 1
                            atlen += 1
                    if _base(dna2, n_rc) == "a":
                        tlen_rc = 0
                        talen_rc = 0
                        if _base(dna2, n_rc - 1) == "t":
                            alen_rc = 0
                            atlen_rc = 0
                        else:
                            alen_rc += 1
                            atlen_rc += 1
                    if _base(dna, n) == "t":
                        if talen < alen:
                            talen += 1
                            atlen += 1
                        else:
                            tlen += 1
                            talen = 0
                            atlen = 0
                            alen = 0
                    if _base(dna2, n_rc) == "t":
                        if talen_rc < alen_rc:
                            talen_rc += 1
                            atlen_rc += 1
                        else:
                            tlen_rc += 1
                            talen_rc = 0
                            atlen_rc = 0
                            alen_rc = 0
                    if max_atlen < atlen:
                        max_atlen = atlen
                        max_atend = n
                    if max_tlen < tlen:
                        max_tlen = tlen
                    if max_atlen_rc < atlen_rc:
                        max_atlen_rc = atlen_rc
                        max_atend_rc = n_rc
                    if max_tlen_rc < tlen_rc:
                        max_tlen_rc = tlen_rc
                if (max_atlen - max_tlen) >= min_at or (max_atlen_rc - max_tlen_rc) >= min_at:
                    if (max_atlen - max_tlen) >= (max_atlen_rc - max_tlen_rc):
                        a_center = (float(max_atend) - ((float(max_atlen) - 1) / 2)) + 1
                    else:
                        a_center = total_bases - (float(max_atend_rc) - ((float(max_atlen_rc) - 1) / 2))
                    patrs.append(PotentialBentDNA(a_center=a_center, strt=strt, end=strt + n_as))
            n_as = 0
        i += 1
    return patrs


def find_apr(context: SequenceContext, min_apr: int, max_apr: int, min_atracts: int) -> list[Repeat]:
    """Find A-phased repeats."""

    processed = _get_atracts(context, min_apr, max_apr)
    repeats: list[Repeat] = []
    tracts = 1

    for i in range(0, len(processed) - (min_atracts + 1)):
        dist_to_next = processed[i + 1].a_center - processed[i].a_center
        if 9.9 <= dist_to_next <= 11.1:
            tracts += 1
        else:
            if tracts >= min_atracts:
                repeats.append(
                    Repeat(
                        start=processed[(i - tracts) + 1].strt,
                        loop=0,
                        num=tracts,
                        strand=0,
                        len=tracts,
                        end=processed[i].end - 1,
                    )
                )
            tracts = 1
    return repeats


def is_subset(repeats: list[Repeat], motif: str, max_loop: int, limit: int, dna: str) -> None:
    """Set the subset flag used by the historical output formats."""

    for rep in repeats:
        rep.special = 0
        if motif == "M":
            n_y = 0
            n_r = 0
            for j in range(rep.start, rep.start + rep.len + 1):
                match _base(dna, j):
                    case "a" | "g":
                        n_r += 1
                    case "c" | "t":
                        n_y += 1
            if n_y == 0:
                ypercent = 0
            elif n_r == 0:
                ypercent = 100
            else:
                ypercent = (n_y // n_r) * 100
            if rep.loop < max_loop and ypercent <= limit:
                rep.special = 1
        elif motif == "D":
            if rep.loop <= max_loop:
                rep.special = 1
        elif motif == "Z":
            if rep.loop >= limit:
                rep.special = 1
        elif motif == "I":
            if rep.len >= limit and rep.loop <= max_loop:
                rep.special = 1
