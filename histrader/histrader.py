#!/usr/bin/env python3
"""
HisTrader -- identify Nucleosome-Free Regions (NFRs) from ChIP-Seq of histone
modifications.

Python port of the original Perl implementation:
    Kirbizakis, Yan & Bailey, 2020 -- https://doi.org/10.1101/2020.03.12.989228
    https://github.com/SvenBaileyLab/HisTrader

"""

from __future__ import annotations
from typing import Optional, List, Dict, Tuple

import argparse
import random
import sys

import numpy as np
import pandas as pd
import pyranges as pr
try:
    import pyBigWig
    _PYBIGWIG = True
except ImportError:
    _PYBIGWIG = False

# An interval is a half-open [start, end) pair of ints, BED-style.
Interval = tuple  # (start, end) as ints

BIGWIG_EXTS = (".bw", ".bigwig", ".bigWig", ".BigWig", ".BW")


def is_bigwig(path: str) -> bool:
    return path.lower().endswith(tuple(e.lower() for e in BIGWIG_EXTS))


# --------------------------------------------------------------------------- #
# Small numeric helpers
# --------------------------------------------------------------------------- #
def round_half_up(x: float) -> int:
    """Round half away from zero, matching the Perl `round` sub's intent
    (its `int($a + $a/abs($a*2))` formula is just a convoluted way to do this)."""
    if x >= 0:
        return int(x + 0.5)
    return -int(-x + 0.5)


# --------------------------------------------------------------------------- #
# I/O
# --------------------------------------------------------------------------- #
def load_genome(path: str) -> dict[str, str]:
    """Read a FASTA file into {seq_name: sequence}. The name is the first
    whitespace-delimited token after '>'."""
    genome: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if line.startswith(">"):
                if name is not None:
                    genome[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        genome[name] = "".join(chunks)
    return genome


def get_fasta(genome: dict[str, str], chrom: str, start: int, end: int) -> str:
    """Sequence for a half-open [start, end) interval. Returns '' if the
    chromosome is absent. Negative starts are clamped (the Perl let a negative
    offset silently index from the end of the string)."""
    seq = genome.get(chrom)
    if seq is None:
        return ""
    return seq[max(start, 0):end]


def load_bedgraph(path: str) -> dict[str, dict[str, np.ndarray]]:
    """Load a bedGraph into per-chromosome arrays, sorted by start.

    Returns {chrom: {"starts", "ends", "signal"}} with numpy arrays. track,
    browser and comment lines are skipped.
    """
    raw: dict[str, list[tuple[int, int, float]]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith(("track", "browser", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            raw.setdefault(f[0], []).append((int(f[1]), int(f[2]), float(f[3])))

    out: dict[str, dict[str, np.ndarray]] = {}
    for chrom, rows in raw.items():
        rows.sort(key=lambda r: r[0])
        out[chrom] = {
            "starts": np.fromiter((r[0] for r in rows), np.int64, len(rows)),
            "ends": np.fromiter((r[1] for r in rows), np.int64, len(rows)),
            "signal": np.fromiter((r[2] for r in rows), np.float64, len(rows)),
        }
    return out


class BigWigSignal:
    """Lazy per-chromosome signal source backed by a bigWig file.

    Instead of loading the entire genome upfront, each chromosome's intervals
    are fetched on first access and cached. For a typical ChIP-seq run where
    peaks cover ~5% of the genome this avoids reading 95% of the file.

    The interface is dict-like ({chrom: {"starts", "ends", "signal"}}) so it
    is a drop-in replacement for the bedGraph dict everywhere in process().
    """

    def __init__(self, path: str) -> None:
        if not _PYBIGWIG:
            raise RuntimeError("pyBigWig is not installed. Run: pip install pyBigWig")
        self._path = path
        self._bw = pyBigWig.open(path)
        self._chroms: dict[str, int] = self._bw.chroms()
        self._cache: dict[str, dict[str, np.ndarray]] = {}

    def __contains__(self, chrom: str) -> bool:
        return chrom in self._chroms

    def __getitem__(self, chrom: str) -> dict[str, np.ndarray]:
        if chrom not in self._cache:
            self._cache[chrom] = self._fetch(chrom)
        return self._cache[chrom]

    def _fetch(self, chrom: str) -> dict[str, np.ndarray]:
        length = self._chroms[chrom]
        intervals = self._bw.intervals(chrom, 0, length)
        if not intervals:
            return {
                "starts": np.empty(0, np.int64),
                "ends":   np.empty(0, np.int64),
                "signal": np.empty(0, np.float64),
            }
        starts = np.fromiter((r[0] for r in intervals), np.int64, len(intervals))
        ends   = np.fromiter((r[1] for r in intervals), np.int64, len(intervals))
        signal = np.fromiter((r[2] for r in intervals), np.float64, len(intervals))
        signal = np.where(np.isnan(signal), 0.0, signal)
        return {"starts": starts, "ends": ends, "signal": signal}

    def close(self) -> None:
        self._bw.close()

    def keys(self):
        return self._chroms.keys()


def load_signal(path: str):
    """Return a signal source for either a bedGraph or bigWig.

    For bedGraph: loads everything upfront into a plain dict (has to -- no index).
    For bigWig:   returns a BigWigSignal that fetches each chromosome lazily on
                  first access, so only chromosomes that actually have peaks get
                  read from disk.
    """
    if is_bigwig(path):
        return BigWigSignal(path)
    return load_bedgraph(path)


def load_peaks(path: str) -> dict[str, list[Interval]]:
    """Load a BED peak file into {chrom: [(start, end), ...]} sorted by start.
    Duplicate intervals on a chromosome are dropped, matching the Perl's
    hash-keyed dedup. Header lines (POSITION / #) are skipped."""
    peaks: dict[str, set[Interval]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[0].upper().startswith("POSITION") or f[0].startswith("#"):
                continue
            if len(f) < 3:
                continue
            peaks.setdefault(f[0], set()).add((int(f[1]), int(f[2])))
    return {chrom: sorted(ivs) for chrom, ivs in peaks.items()}


# --------------------------------------------------------------------------- #
# Signal processing
# --------------------------------------------------------------------------- #
def moving_average_centre(values: np.ndarray, window: int) -> np.ndarray:
    """Centred moving average with a window of `window` bins.

    Mirrors the Perl movingAverageCentre exactly, including its edge behaviour:
    near the array ends the window is truncated and the average is taken over
    only the bins that exist (the Perl padded with zeros for the *sum* but
    incremented the *count* only for real bins, which works out to the same
    truncated-window mean). For window <= 1 this is the identity.
    """
    n = len(values)
    if n == 0 or window <= 1:
        return values.astype(np.float64, copy=True)

    before = (window - 1) // 2          # floor((window-1)/2)
    after = window - 1 - before          # ceil((window-1)/2)

    # Cumulative sum with a leading zero lets us take any window sum in O(1).
    csum = np.concatenate(([0.0], np.cumsum(values, dtype=np.float64)))
    idx = np.arange(n)
    lo = np.maximum(idx - before, 0)
    hi = np.minimum(idx + after, n - 1)
    window_sum = csum[hi + 1] - csum[lo]
    window_len = (hi - lo + 1).astype(np.float64)
    return window_sum / window_len


def diff_smoothed(values: np.ndarray, window: int) -> np.ndarray:
    """Smooth with a centred moving average, then take the first difference.
    Length is len(values) - 1 (empty if values has 0 or 1 elements)."""
    ma = moving_average_centre(values, window)
    return np.diff(ma)


# --------------------------------------------------------------------------- #
# Interval operations
# --------------------------------------------------------------------------- #
def merge_close(intervals: list[Interval], gap: int) -> list[Interval]:
    """Merge intervals (assumed sorted by start) that are within `gap` bp of
    each other. Equivalent to the Perl `merge`.

    Vectorized for large lists: uses numpy to find merge boundaries in one pass,
    then rebuilds the merged list with a single Python loop over groups.
    For very short lists (< 8 intervals) the Python path is faster.
    """
    n = len(intervals)
    if n == 0:
        return []
    if n < 8:
        merged: list[Interval] = [intervals[0]]
        for start, end in intervals[1:]:
            last_start, last_end = merged[-1]
            if last_end + gap >= start:
                merged[-1] = (last_start, max(last_end, end))
            else:
                merged.append((start, end))
        return merged

    starts = np.fromiter((iv[0] for iv in intervals), np.int64, n)
    ends   = np.fromiter((iv[1] for iv in intervals), np.int64, n)

    # A new group starts wherever ends[i-1] + gap < starts[i].
    new_group = np.empty(n, dtype=bool)
    new_group[0] = True
    new_group[1:] = (ends[:-1] + gap) < starts[1:]
    group_ids = np.cumsum(new_group) - 1   # 0-based group index per interval

    n_groups = int(group_ids[-1]) + 1
    out: list[Interval] = []
    for g in range(n_groups):
        mask = group_ids == g
        out.append((int(starts[mask][0]), int(ends[mask].max())))
    return out


def merge_overlaps(intervals: list[Interval]) -> list[Interval]:
    """Merge strictly overlapping/adjacent intervals. Equivalent to the Perl
    `mergeOverlaps`, minus its stray debug print."""
    return merge_close(sorted(intervals), gap=0)


def gaps_between(intervals: list[Interval]) -> list[Interval]:
    """The gaps between consecutive intervals -- i.e. NFRs between nucleosomes.
    Equivalent to the Perl `getNucFreeDiff`."""
    out: list[Interval] = []
    for (a_start, a_end), (b_start, b_end) in zip(intervals, intervals[1:]):
        out.append((a_end, b_start))
    return out


def _to_pr(intervals: list[Interval], chrom: str = "x") -> "pr.PyRanges":
    """Convert a list of (start, end) tuples to a single-chromosome PyRanges."""
    df = pd.DataFrame({
        "Chromosome": chrom,
        "Start": [iv[0] for iv in intervals],
        "End": [iv[1] for iv in intervals],
    })
    return pr.PyRanges(df)


def _from_pr(ranges: "pr.PyRanges") -> list[Interval]:
    """Convert a PyRanges back to a sorted list of (start, end) tuples."""
    if len(ranges) == 0:
        return []
    df = ranges.df.sort_values("Start")
    return list(zip(df["Start"].tolist(), df["End"].tolist()))


def interval_union(a: list[Interval], b: list[Interval]) -> list[Interval]:
    """Union of two interval sets via pyranges set_union + merge.

    Replaces the Perl overlap_union which was O(n^2) pairwise-only and computed
    bounding spans rather than a true union.
    """
    if not a and not b:
        return []
    if not a:
        return _from_pr(_to_pr(b).merge())
    if not b:
        return _from_pr(_to_pr(a).merge())
    return _from_pr(_to_pr(a).set_union(_to_pr(b)))


def interval_intersection(a: list[Interval], b: list[Interval]) -> list[Interval]:
    """Intersection of two interval sets via pyranges set_intersect.

    Replaces the Perl overlap_consensus (O(n^2) pairwise, strict comparators
    that dropped shared-endpoint overlaps).
    """
    if not a or not b:
        return []
    return _from_pr(_to_pr(a).set_intersect(_to_pr(b)))


def filter_by_size(intervals: list[Interval], max_size: int) -> list[Interval]:
    """Drop intervals wider than `max_size`. Equivalent to the Perl `filterBed`."""
    return [iv for iv in intervals if (iv[1] - iv[0]) <= max_size]


# --------------------------------------------------------------------------- #
# Per-peak signal extraction
# --------------------------------------------------------------------------- #
def extract_peak_signal(
    bg: dict[str, np.ndarray], peak: Interval
) -> tuple[np.ndarray, np.ndarray]:
    """Single-bp-resolution signal across a peak.

    Returns (signal, positions) where positions[i] = (i, i+1) and signal[i] is
    the bedGraph value at that base. Combines the Perl highRes expansion with
    getProbeInt's overlap selection.

    Vectorized: each bedGraph interval is clipped to the peak, then its signal
    value is tiled for every base it covers via np.repeat -- no Python loop over
    individual bases.
    """
    p_start, p_end = peak
    starts, ends, signal = bg["starts"], bg["ends"], bg["signal"]

    # bedGraph intervals overlapping [p_start, p_end).
    lo = np.searchsorted(ends, p_start, side="right")
    hi = np.searchsorted(starts, p_end, side="left")
    if hi <= lo:
        return np.empty(0, np.float64), np.empty((0, 2), np.int64)

    clipped_starts = np.maximum(starts[lo:hi], p_start)
    clipped_ends   = np.minimum(ends[lo:hi],   p_end)
    widths = clipped_ends - clipped_starts
    mask = widths > 0
    if not mask.any():
        return np.empty(0, np.float64), np.empty((0, 2), np.int64)

    clipped_starts = clipped_starts[mask]
    clipped_ends   = clipped_ends[mask]
    widths         = widths[mask]
    sig_slice      = signal[lo:hi][mask]

    # Tile each signal value for every base it covers.
    hres_sig = np.repeat(sig_slice, widths)

    # Build base-start positions fully vectorized using cumulative offsets.
    # For each clipped interval [s, e), we want arange(s, e). Rather than
    # calling np.arange N times, we build a single sequence of 0..total_bases-1
    # and shift each segment by its interval's start using cumulative widths.
    total = int(widths.sum())
    base_idx = np.arange(total, dtype=np.int64)          # 0,1,...,total-1
    # For each interval i, subtract the cumulative offset so that the first
    # base of interval i lands at clipped_starts[i].
    cum = np.concatenate(([0], np.cumsum(widths[:-1])))  # start index of each interval
    # repeat each interval's start and its cumulative offset
    interval_starts = np.repeat(clipped_starts, widths)
    interval_cum    = np.repeat(cum, widths)
    offsets = interval_starts + (base_idx - interval_cum)
    hres_pos = np.stack([offsets, offsets + 1], axis=1)
    return hres_sig, hres_pos


def rebin_fixed_step(
    signal: np.ndarray, positions: np.ndarray, step: int
) -> tuple[np.ndarray, list[Interval]]:
    """Re-bin single-bp signal to fixed-width `step` bins.

    Each output bin takes the *max* signal of the bases it covers and spans from
    the first base's start to the last base's end. Equivalent to the Perl
    `increaseStep`.

    Vectorized: the signal array is zero-padded to a multiple of `step`, then
    reshaped into (n_bins, step) so np.max runs over entire rows in one call.
    Bin boundaries are read directly from the positions array without a loop.

    FIXED: the Perl's increaseStep closed a bin one iteration late, so the
    boundary base appeared in two bins' step tracking. Each base belongs to
    exactly one bin here.
    """
    n = len(signal)
    if n == 0:
        return np.empty(0, np.float64), []

    n_bins = (n + step - 1) // step
    pad = n_bins * step - n

    # Pad with -inf so padding never wins the max.
    padded = np.concatenate([signal, np.full(pad, -np.inf)])
    new_sig = padded.reshape(n_bins, step).max(axis=1)

    # Bin i spans positions[i*step] to positions[min((i+1)*step-1, n-1)].
    bin_starts = positions[np.arange(n_bins) * step, 0]
    last_idx   = np.minimum(np.arange(1, n_bins + 1) * step, n) - 1
    bin_ends   = positions[last_idx, 1]

    new_step = list(zip(bin_starts.tolist(), bin_ends.tolist()))
    return new_sig, new_step


# --------------------------------------------------------------------------- #
# NFR calling
# --------------------------------------------------------------------------- #
def call_diff(signal: np.ndarray, step_intervals: list[Interval],
              window: int, merge_dist: int) -> tuple[list[Interval], list[Interval]]:
    """DIFF method. Returns (nucleosomes, nfrs).

    Second-order differencing: bins where the 2nd difference is negative are
    flagged, merged into nucleosomes, and the gaps between them are NFRs. The
    `+ 2` index offset matches the Perl diffOrd_2 (two diffs each drop one bin
    and it indexes the trailing edge).
    """
    first = diff_smoothed(signal, window)
    second = diff_smoothed(first, window)

    neg = np.nonzero(second < 0)[0] + 2
    neg = neg[neg < len(step_intervals)]
    if neg.size == 0:
        return [], []

    flagged = [step_intervals[i] for i in neg]
    nucleosomes = merge_close(flagged, merge_dist)
    nfrs = gaps_between(nucleosomes)
    return nucleosomes, nfrs


def call_ma(signal: np.ndarray, step_intervals: list[Interval],
            fast_window: int, slow_window: int,
            merge_dist: int) -> tuple[list[Interval], list[Interval]]:
    """MA method. Returns (nucleosomes, nfrs).

    Bins where the fast moving average exceeds the slow one are nucleosomal;
    those bins are merged into nucleosomes and the gaps between them are NFRs.
    """
    fast = moving_average_centre(signal, fast_window)
    slow = moving_average_centre(signal, slow_window)

    above = np.nonzero(fast > slow)[0]
    if above.size == 0:
        return [], []

    flagged = [step_intervals[i] for i in above]
    nucleosomes = merge_close(flagged, merge_dist)
    nfrs = gaps_between(nucleosomes)
    return nucleosomes, nfrs


def max_signal_in(interval: Interval, signal: np.ndarray,
                  step_intervals: list[Interval]) -> float | None:
    """Max bin signal whose bin lies fully within `interval`."""
    s, e = interval
    vals = [signal[i] for i, (bs, be) in enumerate(step_intervals)
            if bs >= s and be <= e]
    return max(vals) if vals else None


def min_signal_in(interval: Interval, signal: np.ndarray,
                   step_intervals: list[Interval]) -> float | None:
    """Min bin signal whose bin lies fully within `interval`."""
    s, e = interval
    vals = [signal[i] for i, (bs, be) in enumerate(step_intervals)
            if bs >= s and be <= e]
    return min(vals) if vals else None


def pick_max_valley(nucleosomes: list[Interval], nfrs: list[Interval],
                    signal: np.ndarray, step_intervals: list[Interval],
                    use_differential: bool) -> Interval | None:
    """Pick a single 'best' NFR for a peak. Equivalent to the Perl
    getMaxNucValley.

    use_differential=True : the NFR maximising (flanking nucleosome max signal
                            sum) - (NFR min signal).
    use_differential=False: of the two NFRs flanking the single highest
                            nucleosome, the one with the lower minimum signal.
    Returns None if no NFR can be chosen.
    """
    if not nucleosomes or not nfrs:
        return None

    nuc_max = {nuc: max_signal_in(nuc, signal, step_intervals)
               for nuc in nucleosomes}

    if use_differential:
        best_nfr: Interval | None = None
        best_diff = -np.inf
        for k in range(len(nucleosomes) - 1):
            left, right = nucleosomes[k], nucleosomes[k + 1]
            if nuc_max[left] is None or nuc_max[right] is None:
                continue
            # FIXED: the Perl indexed nfrs[k+1] here, but gaps_between() yields
            # exactly one NFR per adjacent nucleosome pair, so the NFR between
            # nucleosomes k and k+1 is nfrs[k]. The Perl's off-by-one made the
            # differential method inspect the wrong valley and never consider
            # the last pair at all.
            if k >= len(nfrs):
                continue
            nfr = nfrs[k]
            nfr_min = min_signal_in(nfr, signal, step_intervals)
            if nfr_min is None:
                nfr_min = 1e9
            differential = nuc_max[left] + nuc_max[right] - nfr_min
            if differential > best_diff:
                best_diff = differential
                best_nfr = nfr
        return best_nfr

    # Non-differential: find the highest nucleosome, then the lower-signal of
    # its two flanking NFRs.
    highest = max(
        (nuc for nuc in nucleosomes if nuc_max[nuc] is not None),
        key=lambda nuc: nuc_max[nuc],
        default=None,
    )
    if highest is None:
        return None
    k = nucleosomes.index(highest)

    left_sig = right_sig = None
    if k > 0 and (k - 1) < len(nfrs):
        left_sig = min_signal_in(nfrs[k - 1], signal, step_intervals)
        if left_sig is None:
            left_sig = 1e9
    if k < len(nucleosomes) - 1 and k < len(nfrs):
        right_sig = min_signal_in(nfrs[k], signal, step_intervals)
        if right_sig is None:
            right_sig = 1e9

    if left_sig is not None and right_sig is not None:
        return nfrs[k - 1] if left_sig < right_sig else nfrs[k]
    if left_sig is not None:
        return nfrs[k - 1]
    if right_sig is not None:
        return nfrs[k]
    return None


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def write_bed(handle, chrom: str, peak: Interval, intervals: list[Interval]) -> None:
    """Write intervals as BED rows. Columns 4-6 carry the parent peak id, the
    total count of intervals for this peak, and this interval's 1-based index,
    matching the Perl printBed."""
    if not intervals:
        return
    total = len(intervals)
    peak_id = f"{chrom}:{peak[0]}-{peak[1]}"
    for idx, (start, end) in enumerate(intervals, start=1):
        handle.write(f"{chrom}\t{start}\t{end}\t{peak_id}\t{total}\t{idx}\n")


def write_fasta(handle, genome: dict[str, str], chrom: str, peak: Interval,
                intervals: list[Interval], trim: bool, trim_size: int) -> None:
    """Write FASTA records for intervals. With trim=True each sequence is
    `trim_size` bp centred on the interval midpoint. Equivalent to printFasta."""
    if not intervals:
        return
    total = len(intervals)
    for idx, (start, end) in enumerate(intervals, start=1):
        if trim:
            mid = round_half_up((start + end) / 2)
            ext = round_half_up(trim_size / 2)
            s, e = mid - ext, mid + ext
        else:
            s, e = start, end
        seq = get_fasta(genome, chrom, s, e)
        handle.write(f">{chrom}:{s}-{e}_{peak[0]}-{peak[1]}_{total}_{idx}\n{seq}\n")


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def process(args: argparse.Namespace) -> None:
    method = args.method.upper()
    if method not in ("MA", "DIFF", "BOTH"):
        sys.exit(f"ERROR: --method must be MA, DIFF or BOTH (got {args.method!r})")
    if args.maxValley and args.randValley:
        sys.exit("ERROR: --maxValley and --randValley cannot be used together.")

    if args.seed is not None:
        random.seed(args.seed)

    merge_dist = args.step * args.mergeMulti
    # Bins per nucleosome -> fast moving-average window.
    fast_window = round_half_up(args.nucSize / args.step)
    # Slow moving-average window: maMulti nucleosomes wide.
    # FIXED: in the Perl, --maMulti was accidentally wired to mergeMulti, so this
    # could never actually be tuned. It is its own parameter here.
    slow_window = round_half_up((args.nucSize * args.maMulti) / args.step)

    print(f"Merge Distance = {merge_dist}")
    print(f"Bins Per Nucleosome (fast MA window) = {fast_window}")
    print(f"Slow MA window = {slow_window}")
    print(f"pMax = {args.pMax}    Filter = {args.filter}")
    if args.maxValley:
        print(f"maxValley = TRUE    useDifferential = {args.useDifferential}")

    genome: dict[str, str] = {}
    if args.genome:
        genome = load_genome(args.genome)

    bedgraph = load_signal(args.bedGraph)
    peaks = load_peaks(args.peaks)

    nfr_path = f"{args.out}.nfr.bed"
    nuc_path = f"{args.out}.nuc.bed"
    miss_path = f"{args.out}.missing.bed"
    print(f"\nOutput Filenames:\n  NFRs                = {nfr_path}")
    print(f"  Nucleosome-occupied = {nuc_path}")
    print(f"  No NFR detected     = {miss_path}")

    nfr_fh = open(nfr_path, "w")
    nuc_fh = open(nuc_path, "w")
    miss_fh = open(miss_path, "w")
    bg_fh = open(f"{args.out}.bedGraph", "w") if args.outBG else None
    nfr_fa = nuc_fa = None
    if genome:
        nfr_fa = open(f"{args.out}.nfr.fa", "w")
        nuc_fa = open(f"{args.out}.nuc.fa", "w")
        print(f"  NFR sequences       = {args.out}.nfr.fa")
        print(f"  Nucleosome seqs     = {args.out}.nuc.fa")

    try:
        for chrom in sorted(peaks):
            if chrom not in bedgraph:
                continue
            bg = bedgraph[chrom]
            for peak in peaks[chrom]:
                if (peak[1] - peak[0]) < args.minSize:
                    continue

                signal, positions = extract_peak_signal(bg, peak)
                if signal.size == 0:
                    continue

                new_signal, step_intervals = rebin_fixed_step(
                    signal, positions, args.step)
                if new_signal.size == 0:
                    continue

                # Zero out bins below the pMax threshold.
                if args.pMax > 0:
                    cutoff = round_half_up(float(new_signal.max()) * args.pMax)
                    new_signal = np.where(new_signal < cutoff, 0.0, new_signal)

                diff_nuc: list[Interval] = []
                diff_nfr: list[Interval] = []
                ma_nuc: list[Interval] = []
                ma_nfr: list[Interval] = []

                if method in ("DIFF", "BOTH"):
                    diff_nuc, diff_nfr = call_diff(
                        new_signal, step_intervals, fast_window, merge_dist)
                    if method == "DIFF" and not diff_nfr:
                        miss_fh.write(f"{chrom}\t{peak[0]}\t{peak[1]}\n")

                if method in ("MA", "BOTH"):
                    ma_nuc, ma_nfr = call_ma(
                        new_signal, step_intervals,
                        fast_window, slow_window, merge_dist)
                    if method == "MA" and not ma_nfr:
                        miss_fh.write(f"{chrom}\t{peak[0]}\t{peak[1]}\n")

                # Resolve to a single (nucleosomes, nfrs) pair per method.
                if method == "DIFF":
                    nucleosomes, nfrs = diff_nuc, diff_nfr
                elif method == "MA":
                    nucleosomes, nfrs = ma_nuc, ma_nfr
                else:  # BOTH
                    if diff_nfr and ma_nfr:
                        nfrs = interval_intersection(diff_nfr, ma_nfr)
                        nucleosomes = interval_union(diff_nuc, ma_nuc)
                    else:
                        miss_fh.write(f"{chrom}\t{peak[0]}\t{peak[1]}\n")
                        nfrs, nucleosomes = [], []

                if not nfrs:
                    if bg_fh is not None:
                        for (bs, be), v in zip(step_intervals, new_signal):
                            bg_fh.write(f"{chrom}\t{bs}\t{be}\t{v}\n")
                    continue

                nfrs = filter_by_size(nfrs, args.filter)

                if args.maxValley:
                    best = pick_max_valley(
                        nucleosomes, nfrs, new_signal, step_intervals,
                        args.useDifferential)
                    nfrs = [best] if best is not None else []
                elif args.randValley and nfrs:
                    nfrs = [random.choice(nfrs)]

                write_bed(nfr_fh, chrom, peak, nfrs)
                write_bed(nuc_fh, chrom, peak, nucleosomes)
                if genome:
                    write_fasta(nfr_fa, genome, chrom, peak, nfrs,
                                args.trim, args.trimSize)
                    write_fasta(nuc_fa, genome, chrom, peak, nucleosomes,
                                args.trim, args.trimSize)

                if bg_fh is not None:
                    for (bs, be), v in zip(step_intervals, new_signal):
                        bg_fh.write(f"{chrom}\t{bs}\t{be}\t{v}\n")
    finally:
        nfr_fh.close()
        nuc_fh.close()
        miss_fh.close()
        if bg_fh is not None:
            bg_fh.close()
        if nfr_fa is not None:
            nfr_fa.close()
        if nuc_fa is not None:
            nuc_fa.close()
        if isinstance(bedgraph, BigWigSignal):
            bedgraph.close()

    print("\nFINISHED!")


def build_parser() -> argparse.ArgumentParser:
    from histrader import __version__
    class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                         argparse.RawDescriptionHelpFormatter):
        pass

    p = argparse.ArgumentParser(
        prog="histrader.py",
        description=HEADER + "\nIdentify Nucleosome-Free Regions from ChIP-Seq "
                    "of histone modifications.",
        formatter_class=_HelpFormatter,
    )
    p.add_argument("--version", action="version",
                   version=f"HisTrader {__version__}")
    p.add_argument("--signal", required=True, dest="bedGraph", metavar="SIGNAL",
                   help="ChIP-Seq signal file in bedGraph or bigWig format (.bw/.bigwig).")
    p.add_argument("--peaks", required=True,
                   help="Broad peak file in BED format.")
    p.add_argument("--genome",
                   help="Genome FASTA; enables DNA sequence extraction for NFRs.")
    p.add_argument("--trim", action="store_true",
                   help="Trim extracted FASTA sequences (needs --genome, "
                        "--trimSize).")
    p.add_argument("--trimSize", type=int, default=100,
                   help="Length (bp) of trimmed FASTA sequences, centred on "
                        "each NFR.")
    p.add_argument("--out", default="Histrader", help="Output file prefix.")
    p.add_argument("--method", default="BOTH", choices=["MA", "DIFF", "BOTH"],
                   help="NFR detection method.")
    p.add_argument("--step", type=int, default=25,
                   help="Fixed step size (bp) for the converted signal.")
    p.add_argument("--minSize", type=int, default=500,
                   help="Minimum peak size (bp) for NFR calling.")
    p.add_argument("--nucSize", type=int, default=150,
                   help="Estimated nucleosome size (bp); should be divisible "
                        "by --step.")
    p.add_argument("--mergeMulti", type=int, default=3,
                   help="Step multiplier for merging (merge distance = "
                        "mergeMulti * step).")
    p.add_argument("--maMulti", type=int, default=3,
                   help="Nucleosome multiplier for the slow moving average.")
    p.add_argument("--pMax", type=float, default=0.0,
                   help="Fraction of peak max signal used as a zeroing "
                        "threshold (0 = off).")
    p.add_argument("--filter", type=int, default=1000,
                   help="Discard NFRs wider than this (bp).")
    p.add_argument("--outBG", action="store_true",
                   help="Also output the fixed-step signal within peaks "
                        "(bedGraph).")
    p.add_argument("--maxValley", action="store_true",
                   help="Keep only the NFR at the max peak region per peak.")
    p.add_argument("--useDifferential", action="store_true",
                   help="With --maxValley, use the differential-based "
                        "selection.")
    p.add_argument("--randValley", action="store_true",
                   help="Keep only one random NFR per peak.")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed, for reproducible --randValley output.")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress the banner header.")
    return p


HEADER = """
########################################################################################################
##                                                                                                    ##
##    HISTRADER: A tool to identify nucleosome free regions from ChIP-Seq of Histone Modifications    ##
##                                                                                                    ##
##                  Written by Eftyhios Kirbizakis, Yifei Yan, and Swneke D. Bailey                   ##
##                               Copyright 2020 Swneke D. Bailey                                      ##
##                                                                                                    ##
########################################################################################################
"""


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if not args.quiet:
        print(HEADER)
    print(f"Identifying valleys in {args.bedGraph} at positions in {args.peaks}\n")
    process(args)


if __name__ == "__main__":
    main()