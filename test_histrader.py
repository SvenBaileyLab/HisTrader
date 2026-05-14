"""
pytest test suite for histrader.py

Tests are grouped into sections matching the module's structure:
  - Numeric helpers
  - I/O (load_genome, get_fasta, load_bedgraph, load_peaks)
  - Signal processing (moving_average_centre, diff_smoothed)
  - Interval operations (merge_close, gaps_between, union, intersection, filter)
  - Per-peak signal extraction & rebinning
  - NFR calling (call_diff, call_ma, pick_max_valley)
  - Integration (process() end-to-end via tmp files)
"""

import math
import os
import random
import textwrap
from pathlib import Path

import numpy as np
import pytest

import histrader as H


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────
@pytest.fixture()
def tmp(tmp_path):
    """Just expose tmp_path with a short name."""
    return tmp_path


@pytest.fixture()
def two_hump_bg(tmp):
    """bedGraph for chr1 0..3000 with two gaussian humps at 800 and 2200."""
    path = tmp / "two_hump.bedGraph"
    with path.open("w") as f:
        for start in range(0, 3000, 10):
            x = start + 5
            sig = round(
                100 * math.exp(-((x - 800) ** 2) / (2 * 250 ** 2))
                + 100 * math.exp(-((x - 2200) ** 2) / (2 * 250 ** 2))
                + 2,
                2,
            )
            f.write(f"chr1\t{start}\t{start+10}\t{sig}\n")
    return str(path)


@pytest.fixture()
def two_hump_peaks(tmp):
    path = tmp / "peaks.bed"
    path.write_text("chr1\t100\t2900\n")
    return str(path)


@pytest.fixture()
def simple_genome(tmp):
    path = tmp / "genome.fa"
    path.write_text(">chr1\n" + "ACGT" * 1000 + "\n")
    return str(path)


# ─────────────────────────────────────────────────────────────────────────────
# Numeric helpers
# ─────────────────────────────────────────────────────────────────────────────
class TestRoundHalfUp:
    def test_positive_half(self):
        assert H.round_half_up(0.5) == 1

    def test_negative_half(self):
        assert H.round_half_up(-0.5) == -1

    def test_zero(self):
        assert H.round_half_up(0.0) == 0

    def test_integer(self):
        assert H.round_half_up(3.0) == 3

    def test_just_below_half(self):
        assert H.round_half_up(2.4999) == 2

    def test_just_above_half(self):
        assert H.round_half_up(2.5001) == 3


# ─────────────────────────────────────────────────────────────────────────────
# I/O
# ─────────────────────────────────────────────────────────────────────────────
class TestLoadGenome:
    def test_basic(self, tmp):
        p = tmp / "g.fa"
        p.write_text(">chr1\nACGT\nACGT\n>chr2\nTTTT\n")
        g = H.load_genome(str(p))
        assert g["chr1"] == "ACGTACGT"
        assert g["chr2"] == "TTTT"

    def test_header_with_extra_tokens(self, tmp):
        p = tmp / "g.fa"
        p.write_text(">chr1 extra info here\nAAAA\n")
        g = H.load_genome(str(p))
        assert "chr1" in g

    def test_blank_and_comment_lines_skipped(self, tmp):
        p = tmp / "g.fa"
        p.write_text(">chr1\n# comment\n\nACGT\n")
        g = H.load_genome(str(p))
        assert g["chr1"] == "ACGT"


class TestGetFasta:
    def test_basic_slice(self):
        g = {"chr1": "ACGTACGT"}
        assert H.get_fasta(g, "chr1", 0, 4) == "ACGT"

    def test_mid_slice(self):
        g = {"chr1": "ACGTACGT"}
        assert H.get_fasta(g, "chr1", 2, 6) == "GTAC"

    def test_missing_chrom(self):
        assert H.get_fasta({}, "chrX", 0, 10) == ""

    def test_negative_start_clamped(self):
        # Perl bug: negative start would silently index from end of string.
        g = {"chr1": "ACGTACGT"}
        seq = H.get_fasta(g, "chr1", -5, 4)
        assert seq == "ACGT"  # clamped to 0


class TestLoadBedgraph:
    def test_basic(self, tmp):
        p = tmp / "a.bedGraph"
        p.write_text("chr1\t0\t10\t5.0\nchr1\t10\t20\t3.0\n")
        bg = H.load_bedgraph(str(p))
        assert "chr1" in bg
        np.testing.assert_array_equal(bg["chr1"]["starts"], [0, 10])
        np.testing.assert_array_equal(bg["chr1"]["signal"], [5.0, 3.0])

    def test_track_header_skipped(self, tmp):
        p = tmp / "b.bedGraph"
        p.write_text("track type=bedGraph\nchr1\t0\t10\t1.0\n")
        bg = H.load_bedgraph(str(p))
        assert len(bg["chr1"]["starts"]) == 1

    def test_multi_chrom(self, tmp):
        p = tmp / "c.bedGraph"
        p.write_text("chr1\t0\t10\t1.0\nchr2\t0\t10\t2.0\n")
        bg = H.load_bedgraph(str(p))
        assert set(bg.keys()) == {"chr1", "chr2"}

    def test_sorted_by_start(self, tmp):
        p = tmp / "d.bedGraph"
        p.write_text("chr1\t20\t30\t1.0\nchr1\t0\t10\t2.0\n")
        bg = H.load_bedgraph(str(p))
        assert list(bg["chr1"]["starts"]) == [0, 20]


class TestLoadPeaks:
    def test_basic(self, tmp):
        p = tmp / "peaks.bed"
        p.write_text("chr1\t100\t500\nchr2\t0\t300\n")
        peaks = H.load_peaks(str(p))
        assert (100, 500) in peaks["chr1"]
        assert (0, 300) in peaks["chr2"]

    def test_duplicates_dropped(self, tmp):
        p = tmp / "peaks.bed"
        p.write_text("chr1\t100\t500\nchr1\t100\t500\n")
        peaks = H.load_peaks(str(p))
        assert len(peaks["chr1"]) == 1

    def test_sorted(self, tmp):
        p = tmp / "peaks.bed"
        p.write_text("chr1\t500\t800\nchr1\t100\t400\n")
        peaks = H.load_peaks(str(p))
        assert peaks["chr1"][0][0] < peaks["chr1"][1][0]

    def test_header_skipped(self, tmp):
        p = tmp / "peaks.bed"
        p.write_text("POSITION\t0\t0\nchr1\t100\t500\n")
        peaks = H.load_peaks(str(p))
        assert "POSITION" not in peaks


# ─────────────────────────────────────────────────────────────────────────────
# Signal processing
# ─────────────────────────────────────────────────────────────────────────────
class TestMovingAverageCentre:
    def test_flat_signal_unchanged(self):
        sig = np.ones(10)
        result = H.moving_average_centre(sig, 3)
        np.testing.assert_allclose(result, np.ones(10))

    def test_window_1_identity(self):
        sig = np.array([1.0, 2.0, 3.0, 4.0])
        np.testing.assert_allclose(H.moving_average_centre(sig, 1), sig)

    def test_known_interior(self):
        # window=3: interior value is mean of 3 neighbours
        sig = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = H.moving_average_centre(sig, 3)
        assert result[2] == pytest.approx(3.0)

    def test_edge_truncated_window(self):
        # At index 0 with window=3, only 2 values available (index 0,1)
        sig = np.array([4.0, 2.0, 2.0])
        result = H.moving_average_centre(sig, 3)
        assert result[0] == pytest.approx((4.0 + 2.0) / 2)

    def test_empty_input(self):
        result = H.moving_average_centre(np.array([]), 3)
        assert len(result) == 0

    def test_output_same_length(self):
        sig = np.random.rand(50)
        assert len(H.moving_average_centre(sig, 7)) == 50

    def test_smoothing_reduces_variance(self):
        rng = np.random.default_rng(0)
        sig = rng.standard_normal(200)
        smoothed = H.moving_average_centre(sig, 15)
        assert smoothed.std() < sig.std()


class TestDiffSmoothed:
    def test_length_is_n_minus_1(self):
        sig = np.ones(10)
        assert len(H.diff_smoothed(sig, 3)) == 9

    def test_constant_signal_gives_zero_diff(self):
        sig = np.ones(20)
        result = H.diff_smoothed(sig, 5)
        np.testing.assert_allclose(result, 0, atol=1e-10)

    def test_monotone_increasing_positive_diff(self):
        sig = np.arange(20, dtype=float)
        result = H.diff_smoothed(sig, 1)
        assert (result > 0).all()


# ─────────────────────────────────────────────────────────────────────────────
# Interval operations
# ─────────────────────────────────────────────────────────────────────────────
class TestMergeClose:
    def test_non_overlapping_unchanged(self):
        ivs = [(0, 10), (20, 30)]
        assert H.merge_close(ivs, gap=5) == [(0, 10), (20, 30)]

    def test_overlapping_merged(self):
        ivs = [(0, 15), (10, 25)]
        assert H.merge_close(ivs, gap=0) == [(0, 25)]

    def test_within_gap_merged(self):
        ivs = [(0, 10), (15, 25)]  # 5 bp gap
        assert H.merge_close(ivs, gap=5) == [(0, 25)]

    def test_empty(self):
        assert H.merge_close([], gap=0) == []

    def test_single(self):
        assert H.merge_close([(5, 10)], gap=0) == [(5, 10)]

    def test_chain_of_three(self):
        ivs = [(0, 10), (10, 20), (20, 30)]
        assert H.merge_close(ivs, gap=0) == [(0, 30)]


class TestGapsBetween:
    def test_single_gap(self):
        nucs = [(100, 300), (500, 700)]
        assert H.gaps_between(nucs) == [(300, 500)]

    def test_two_gaps(self):
        nucs = [(0, 100), (200, 300), (400, 500)]
        assert H.gaps_between(nucs) == [(100, 200), (300, 400)]

    def test_single_interval_no_gap(self):
        assert H.gaps_between([(0, 100)]) == []

    def test_empty(self):
        assert H.gaps_between([]) == []


class TestIntervalUnion:
    def test_non_overlapping(self):
        a = [(0, 10)]
        b = [(20, 30)]
        result = H.interval_union(a, b)
        assert set(result) == {(0, 10), (20, 30)}

    def test_overlapping_merged(self):
        a = [(0, 20)]
        b = [(10, 30)]
        assert H.interval_union(a, b) == [(0, 30)]

    def test_chain_of_three_collapses(self):
        # The Perl overlap_union would NOT fully collapse this.
        a = [(0, 15), (25, 40)]
        b = [(10, 30)]
        result = H.interval_union(a, b)
        assert result == [(0, 40)]

    def test_empty_a(self):
        assert H.interval_union([], [(5, 10)]) == [(5, 10)]

    def test_empty_b(self):
        assert H.interval_union([(5, 10)], []) == [(5, 10)]

    def test_both_empty(self):
        assert H.interval_union([], []) == []


class TestIntervalIntersection:
    def test_no_overlap(self):
        assert H.interval_intersection([(0, 10)], [(20, 30)]) == []

    def test_partial_overlap(self):
        assert H.interval_intersection([(0, 20)], [(10, 30)]) == [(10, 20)]

    def test_containment(self):
        assert H.interval_intersection([(0, 100)], [(20, 40)]) == [(20, 40)]

    def test_shared_endpoint_not_included(self):
        # Half-open: [0,10) and [10,20) share a point but don't overlap.
        assert H.interval_intersection([(0, 10)], [(10, 20)]) == []

    def test_empty_a(self):
        assert H.interval_intersection([], [(0, 10)]) == []

    def test_multi_interval(self):
        a = [(0, 20), (40, 60)]
        b = [(10, 50)]
        result = H.interval_intersection(a, b)
        assert (10, 20) in result
        assert (40, 50) in result


class TestFilterBySize:
    def test_keeps_small(self):
        ivs = [(0, 100), (200, 250)]
        assert H.filter_by_size(ivs, 1000) == ivs

    def test_drops_large(self):
        ivs = [(0, 2000), (3000, 3500)]
        assert H.filter_by_size(ivs, 1000) == [(3000, 3500)]

    def test_boundary_inclusive(self):
        # Exactly max_size should be kept.
        assert H.filter_by_size([(0, 100)], 100) == [(0, 100)]

    def test_empty(self):
        assert H.filter_by_size([], 500) == []


# ─────────────────────────────────────────────────────────────────────────────
# Per-peak signal extraction & rebinning
# ─────────────────────────────────────────────────────────────────────────────
class TestExtractPeakSignal:
    def _bg(self, rows):
        """Build a bg dict from [(start, end, signal)] rows (already sorted)."""
        rows = sorted(rows)
        return {
            "starts": np.array([r[0] for r in rows], np.int64),
            "ends":   np.array([r[1] for r in rows], np.int64),
            "signal": np.array([r[2] for r in rows], np.float64),
        }

    def test_single_interval(self):
        bg = self._bg([(0, 5, 3.0)])
        sig, pos = H.extract_peak_signal(bg, (0, 5))
        assert len(sig) == 5
        np.testing.assert_array_equal(sig, [3.0] * 5)

    def test_clips_to_peak(self):
        bg = self._bg([(0, 100, 7.0)])
        sig, pos = H.extract_peak_signal(bg, (20, 30))
        assert len(sig) == 10
        assert pos[0, 0] == 20
        assert pos[-1, 1] == 30

    def test_no_overlap_returns_empty(self):
        bg = self._bg([(0, 10, 1.0)])
        sig, pos = H.extract_peak_signal(bg, (50, 100))
        assert sig.size == 0

    def test_positions_are_contiguous_bp(self):
        bg = self._bg([(0, 10, 1.0)])
        sig, pos = H.extract_peak_signal(bg, (0, 10))
        # Each position row should be (k, k+1).
        np.testing.assert_array_equal(pos[:, 0], np.arange(10))
        np.testing.assert_array_equal(pos[:, 1], np.arange(1, 11))

    def test_multiple_intervals(self):
        bg = self._bg([(0, 5, 1.0), (5, 10, 2.0)])
        sig, _ = H.extract_peak_signal(bg, (0, 10))
        assert list(sig[:5]) == [1.0] * 5
        assert list(sig[5:]) == [2.0] * 5


class TestRebinFixedStep:
    def _pos(self, n):
        """Trivial 1-bp position array for n bases starting at 0."""
        starts = np.arange(n, dtype=np.int64)
        return np.stack([starts, starts + 1], axis=1)

    def test_single_bin(self):
        sig = np.array([1.0, 3.0, 2.0])
        pos = self._pos(3)
        new_sig, steps = H.rebin_fixed_step(sig, pos, 10)
        assert len(new_sig) == 1
        assert new_sig[0] == 3.0

    def test_exact_multiple(self):
        sig = np.array([1.0, 3.0, 2.0, 4.0])
        pos = self._pos(4)
        new_sig, steps = H.rebin_fixed_step(sig, pos, 2)
        assert list(new_sig) == [3.0, 4.0]
        assert steps[0] == (0, 2)
        assert steps[1] == (2, 4)

    def test_partial_last_bin(self):
        sig = np.array([1.0, 2.0, 5.0])
        pos = self._pos(3)
        new_sig, steps = H.rebin_fixed_step(sig, pos, 2)
        assert len(new_sig) == 2
        assert new_sig[1] == 5.0   # last bin is just one element

    def test_empty(self):
        sig = np.empty(0)
        pos = np.empty((0, 2), np.int64)
        ns, st = H.rebin_fixed_step(sig, pos, 25)
        assert ns.size == 0
        assert st == []

    def test_step_1_identity(self):
        sig = np.array([1.0, 2.0, 3.0])
        pos = self._pos(3)
        ns, _ = H.rebin_fixed_step(sig, pos, 1)
        np.testing.assert_array_equal(ns, sig)


# ─────────────────────────────────────────────────────────────────────────────
# NFR calling
# ─────────────────────────────────────────────────────────────────────────────
class TestCallDiff:
    """DIFF method on synthetic two-hump signal."""

    def _two_hump_signal(self, step=25):
        n = 3000 // step
        xs = (np.arange(n) + 0.5) * step
        sig = (
            100 * np.exp(-((xs - 800) ** 2) / (2 * 250 ** 2))
            + 100 * np.exp(-((xs - 2200) ** 2) / (2 * 250 ** 2))
        )
        steps = [(i * step, (i + 1) * step) for i in range(n)]
        return sig, steps

    def test_detects_two_nucleosomes(self):
        sig, steps = self._two_hump_signal()
        nucs, nfrs = H.call_diff(sig, steps, window=6, merge_dist=75)
        assert len(nucs) == 2

    def test_detects_one_nfr_between_humps(self):
        sig, steps = self._two_hump_signal()
        nucs, nfrs = H.call_diff(sig, steps, window=6, merge_dist=75)
        assert len(nfrs) == 1

    def test_nfr_between_nucleosomes(self):
        sig, steps = self._two_hump_signal()
        nucs, nfrs = H.call_diff(sig, steps, window=6, merge_dist=75)
        if nfrs and nucs:
            nfr_centre = (nfrs[0][0] + nfrs[0][1]) / 2
            assert 1200 < nfr_centre < 1800

    def test_flat_signal_no_nfrs(self):
        sig = np.ones(100)
        steps = [(i * 25, (i + 1) * 25) for i in range(100)]
        nucs, nfrs = H.call_diff(sig, steps, window=6, merge_dist=75)
        assert nfrs == []


class TestCallMA:
    """MA method on synthetic two-hump signal."""

    def _two_hump_signal(self, step=25):
        n = 3000 // step
        xs = (np.arange(n) + 0.5) * step
        sig = (
            100 * np.exp(-((xs - 800) ** 2) / (2 * 250 ** 2))
            + 100 * np.exp(-((xs - 2200) ** 2) / (2 * 250 ** 2))
        )
        steps = [(i * step, (i + 1) * step) for i in range(n)]
        return sig, steps

    def test_detects_nucleosomes(self):
        sig, steps = self._two_hump_signal()
        nucs, nfrs = H.call_ma(sig, steps, fast_window=6, slow_window=18,
                                merge_dist=75)
        assert len(nucs) >= 1

    def test_nfr_in_valley(self):
        sig, steps = self._two_hump_signal()
        nucs, nfrs = H.call_ma(sig, steps, fast_window=6, slow_window=18,
                                merge_dist=75)
        if nfrs:
            nfr_centre = (nfrs[0][0] + nfrs[0][1]) / 2
            assert 1000 < nfr_centre < 2000

    def test_flat_signal_no_nfrs(self):
        sig = np.ones(100)
        steps = [(i * 25, (i + 1) * 25) for i in range(100)]
        nucs, nfrs = H.call_ma(sig, steps, fast_window=6, slow_window=18,
                                merge_dist=75)
        assert nfrs == []


class TestPickMaxValley:
    def _setup(self):
        """Two nucleosomes [0,100) and [200,300), NFR [100,200)."""
        nucs = [(0, 100), (200, 300)]
        nfrs = [(100, 200)]
        steps = [(i, i + 1) for i in range(300)]
        # Signal: high at nuc centres, low in NFR.
        sig = np.zeros(300)
        sig[0:100] = 80.0
        sig[100:200] = 5.0
        sig[200:300] = 70.0
        return nucs, nfrs, np.array(sig), steps

    def test_returns_interval(self):
        nucs, nfrs, sig, steps = self._setup()
        result = H.pick_max_valley(nucs, nfrs, sig, steps, use_differential=False)
        assert result is not None
        assert isinstance(result, tuple) and len(result) == 2

    def test_differential_returns_correct_nfr(self):
        nucs, nfrs, sig, steps = self._setup()
        result = H.pick_max_valley(nucs, nfrs, sig, steps, use_differential=True)
        assert result == (100, 200)

    def test_non_differential_returns_correct_nfr(self):
        nucs, nfrs, sig, steps = self._setup()
        result = H.pick_max_valley(nucs, nfrs, sig, steps, use_differential=False)
        assert result == (100, 200)

    def test_no_nucleosomes_returns_none(self):
        result = H.pick_max_valley([], [(10, 20)], np.ones(30),
                                   [(i, i+1) for i in range(30)], False)
        assert result is None

    def test_no_nfrs_returns_none(self):
        result = H.pick_max_valley([(0, 10)], [], np.ones(10),
                                   [(i, i+1) for i in range(10)], False)
        assert result is None


# ─────────────────────────────────────────────────────────────────────────────
# Integration tests
# ─────────────────────────────────────────────────────────────────────────────
class TestProcessEndToEnd:

    def _args(self, tmp, extra=None, method="BOTH"):
        """Return a Namespace with required args pointing at tmp test data."""
        import argparse
        ns = argparse.Namespace(
            bedGraph=str(tmp / "signal.bedGraph"),
            peaks=str(tmp / "peaks.bed"),
            genome=None,
            trim=False,
            trimSize=100,
            out=str(tmp / "out"),
            method=method,
            step=25,
            minSize=500,
            nucSize=150,
            mergeMulti=3,
            maMulti=3,
            pMax=0.0,
            filter=1000,
            outBG=False,
            maxValley=False,
            useDifferential=False,
            randValley=False,
            seed=None,
        )
        if extra:
            for k, v in extra.items():
                setattr(ns, k, v)
        return ns

    def _write_two_hump(self, tmp):
        import math
        bg = tmp / "signal.bedGraph"
        with bg.open("w") as f:
            for start in range(0, 3000, 10):
                x = start + 5
                sig = round(
                    100 * math.exp(-((x - 800) ** 2) / (2 * 250 ** 2))
                    + 100 * math.exp(-((x - 2200) ** 2) / (2 * 250 ** 2))
                    + 2, 2
                )
                f.write(f"chr1\t{start}\t{start+10}\t{sig}\n")
        (tmp / "peaks.bed").write_text("chr1\t100\t2900\n")

    def test_diff_produces_nfr(self, tmp):
        self._write_two_hump(tmp)
        H.process(self._args(tmp, method="DIFF"))
        nfr = (tmp / "out.nfr.bed").read_text().strip()
        assert nfr != ""
        fields = nfr.split("\t")
        assert fields[0] == "chr1"

    def test_ma_produces_nfr(self, tmp):
        self._write_two_hump(tmp)
        H.process(self._args(tmp, method="MA"))
        assert (tmp / "out.nfr.bed").stat().st_size > 0

    def test_both_produces_nfr(self, tmp):
        self._write_two_hump(tmp)
        H.process(self._args(tmp, method="BOTH"))
        assert (tmp / "out.nfr.bed").stat().st_size > 0

    def test_output_files_created(self, tmp):
        self._write_two_hump(tmp)
        H.process(self._args(tmp))
        for suffix in (".nfr.bed", ".nuc.bed", ".missing.bed"):
            assert (tmp / f"out{suffix}").exists()

    def test_outbg_file_created(self, tmp):
        self._write_two_hump(tmp)
        H.process(self._args(tmp, extra={"outBG": True}))
        assert (tmp / "out.bedGraph").exists()
        assert (tmp / "out.bedGraph").stat().st_size > 0

    def test_genome_fasta_output(self, tmp):
        self._write_two_hump(tmp)
        gf = tmp / "genome.fa"
        gf.write_text(">chr1\n" + "ACGT" * 1000 + "\n")
        H.process(self._args(tmp, extra={"genome": str(gf)}))
        nfr_fa = (tmp / "out.nfr.fa").read_text()
        assert nfr_fa.startswith(">")

    def test_peak_too_small_skipped(self, tmp):
        import math
        bg = tmp / "signal.bedGraph"
        with bg.open("w") as f:
            for start in range(0, 400, 10):
                f.write(f"chr1\t{start}\t{start+10}\t5.0\n")
        (tmp / "peaks.bed").write_text("chr1\t0\t300\n")  # < minSize=500
        H.process(self._args(tmp))
        assert (tmp / "out.nfr.bed").read_text() == ""

    def test_maxvalley_single_nfr_per_peak(self, tmp):
        self._write_two_hump(tmp)
        H.process(self._args(tmp, method="DIFF",
                              extra={"maxValley": True}))
        lines = (tmp / "out.nfr.bed").read_text().strip().splitlines()
        assert len(lines) == 1

    def test_randvalley_single_nfr_per_peak(self, tmp):
        self._write_two_hump(tmp)
        H.process(self._args(tmp, method="DIFF",
                              extra={"randValley": True, "seed": 0}))
        lines = (tmp / "out.nfr.bed").read_text().strip().splitlines()
        assert len(lines) == 1

    def test_randvalley_reproducible_with_seed(self, tmp):
        self._write_two_hump(tmp)
        H.process(self._args(tmp, method="DIFF",
                              extra={"randValley": True, "seed": 99,
                                     "out": str(tmp / "run1")}))
        H.process(self._args(tmp, method="DIFF",
                              extra={"randValley": True, "seed": 99,
                                     "out": str(tmp / "run2")}))
        assert (tmp / "run1.nfr.bed").read_text() == \
               (tmp / "run2.nfr.bed").read_text()

    def test_conflicting_flags_raises(self, tmp):
        self._write_two_hump(tmp)
        with pytest.raises(SystemExit):
            H.process(self._args(tmp, extra={"maxValley": True,
                                              "randValley": True}))

    def test_chrom_absent_in_bedgraph_skipped(self, tmp):
        """Peaks on a chromosome with no bedGraph data produce no output."""
        bg = tmp / "signal.bedGraph"
        bg.write_text("chr1\t0\t3000\t5.0\n")
        (tmp / "peaks.bed").write_text("chr9\t0\t2000\n")
        H.process(self._args(tmp))
        assert (tmp / "out.nfr.bed").read_text() == ""

    def test_nuc_bed_format(self, tmp):
        """NUC bed rows: chrom, start, end, peak_id, total, index."""
        self._write_two_hump(tmp)
        H.process(self._args(tmp, method="DIFF"))
        lines = (tmp / "out.nuc.bed").read_text().strip().splitlines()
        assert len(lines) >= 1
        for line in lines:
            fields = line.split("\t")
            assert len(fields) == 6
            assert int(fields[1]) < int(fields[2])   # start < end

    def test_pmax_zeros_low_signal(self, tmp):
        """pMax > 0 should zero bins below the threshold and may suppress NFRs
        from genuinely noisy / low-signal peaks."""
        import math
        bg = tmp / "signal.bedGraph"
        with bg.open("w") as f:
            # Very flat signal -- all near 5.0.
            for start in range(0, 3000, 10):
                f.write(f"chr1\t{start}\t{start+10}\t5.0\n")
        (tmp / "peaks.bed").write_text("chr1\t100\t2900\n")
        H.process(self._args(tmp, method="DIFF",
                              extra={"pMax": 0.99}))  # nearly everything zeroed
        # We don't assert a specific output, just that it runs cleanly.
        assert (tmp / "out.nfr.bed").exists()


# ─────────────────────────────────────────────────────────────────────────────
# bigWig support
# ─────────────────────────────────────────────────────────────────────────────
class TestIsBigwig:
    def test_bw_extension(self):
        assert H.is_bigwig("signal.bw")

    def test_bigwig_extension(self):
        assert H.is_bigwig("signal.bigwig")

    def test_mixed_case(self):
        assert H.is_bigwig("signal.BigWig")

    def test_bedgraph_not_bigwig(self):
        assert not H.is_bigwig("signal.bedGraph")

    def test_bed_not_bigwig(self):
        assert not H.is_bigwig("peaks.bed")


class TestLoadBigwig:
    def _make_bw(self, tmp, rows):
        """Write a minimal bigWig from [(chrom, start, end, value)] rows."""
        import pyBigWig
        sizes = {}
        for c, s, e, v in rows:
            sizes[c] = max(sizes.get(c, 0), e + 1)
        path = str(tmp / "test.bw")
        bw = pyBigWig.open(path, "w")
        bw.addHeader(list(sizes.items()))
        for c, s, e, v in sorted(rows):
            bw.addEntries([c], [s], ends=[e], values=[float(v)])
        bw.close()
        return path

    def test_basic_load(self, tmp):
        path = self._make_bw(tmp, [("chr1", 0, 10, 5.0), ("chr1", 10, 20, 3.0)])
        bg = H.load_signal(path)
        assert "chr1" in bg
        assert len(bg["chr1"]["starts"]) == 2

    def test_values_match(self, tmp):
        path = self._make_bw(tmp, [("chr1", 0, 10, 7.0)])
        bg = H.load_signal(path)
        assert float(bg["chr1"]["signal"][0]) == pytest.approx(7.0, rel=1e-4)

    def test_multi_chrom(self, tmp):
        path = self._make_bw(tmp, [("chr1", 0, 10, 1.0), ("chr2", 0, 10, 2.0)])
        bg = H.load_signal(path)
        assert set(bg.keys()) == {"chr1", "chr2"}

    def test_sorted_by_start(self, tmp):
        path = self._make_bw(tmp, [("chr1", 20, 30, 1.0), ("chr1", 0, 10, 2.0)])
        bg = H.load_signal(path)
        assert list(bg["chr1"]["starts"]) == [0, 20]

    def test_nan_replaced_with_zero(self, tmp):
        """Gaps in bigWig coverage return NaN from pyBigWig; we replace with 0."""
        import pyBigWig
        path = str(tmp / "gap.bw")
        bw = pyBigWig.open(path, "w")
        bw.addHeader([("chr1", 1000)])
        bw.addEntries(["chr1"], [100], ends=[200], values=[5.0])
        bw.close()
        bg = H.load_signal(path)
        assert not any(np.isnan(bg["chr1"]["signal"]))

    def test_lazy_only_fetches_accessed_chrom(self, tmp):
        """BigWigSignal should not fetch chr2 if we only access chr1."""
        path = self._make_bw(tmp, [("chr1", 0, 10, 1.0), ("chr2", 0, 10, 2.0)])
        bws = H.load_signal(path)
        assert isinstance(bws, H.BigWigSignal)
        assert "chr2" not in bws._cache   # not yet fetched
        _ = bws["chr1"]
        assert "chr1" in bws._cache
        assert "chr2" not in bws._cache   # still not fetched
        bws.close()

    def test_lazy_caches_on_second_access(self, tmp):
        path = self._make_bw(tmp, [("chr1", 0, 10, 5.0)])
        bws = H.load_signal(path)
        a = bws["chr1"]
        b = bws["chr1"]
        assert a is b   # same object returned from cache
        bws.close()


class TestLoadSignal:
    def test_detects_bedgraph(self, tmp):
        p = tmp / "a.bedGraph"
        p.write_text("chr1\t0\t10\t5.0\n")
        bg = H.load_signal(str(p))
        assert "chr1" in bg

    def test_detects_bigwig(self, tmp):
        import pyBigWig
        path = str(tmp / "a.bw")
        bw = pyBigWig.open(path, "w")
        bw.addHeader([("chr1", 100)])
        bw.addEntries(["chr1"], [0], ends=[50], values=[3.0])
        bw.close()
        bg = H.load_signal(path)
        assert "chr1" in bg


class TestBigwigEndToEnd:
    """Integration: bigWig input should give identical results to bedGraph input
    built from the same data."""

    def _write_bedgraph(self, tmp, rows):
        p = tmp / "signal.bedGraph"
        with p.open("w") as f:
            for c, s, e, v in rows:
                f.write(f"{c}\t{s}\t{e}\t{v}\n")
        return str(p)

    def _write_bigwig(self, tmp, rows):
        import pyBigWig
        sizes = {}
        for c, s, e, v in rows:
            sizes[c] = max(sizes.get(c, 0), e + 1)
        path = str(tmp / "signal.bw")
        bw = pyBigWig.open(path, "w")
        bw.addHeader(list(sizes.items()))
        for c, s, e, v in sorted(rows):
            bw.addEntries([c], [s], ends=[e], values=[float(v)])
        bw.close()
        return path

    def _two_hump_rows(self):
        import math
        rows = []
        for start in range(0, 3000, 10):
            x = start + 5
            sig = round(
                100 * math.exp(-((x - 800) ** 2) / (2 * 250 ** 2))
                + 100 * math.exp(-((x - 2200) ** 2) / (2 * 250 ** 2))
                + 2, 4
            )
            rows.append(("chr1", start, start + 10, sig))
        return rows

    def _args(self, tmp, signal_path, method="BOTH"):
        import argparse
        return argparse.Namespace(
            bedGraph=signal_path,
            peaks=str(tmp / "peaks.bed"),
            genome=None, trim=False, trimSize=100,
            out=str(tmp / "out"),
            method=method, step=25, minSize=500, nucSize=150,
            mergeMulti=3, maMulti=3, pMax=0.0, filter=1000,
            outBG=False, maxValley=False, useDifferential=False,
            randValley=False, seed=None,
        )

    def test_bigwig_matches_bedgraph(self, tmp):
        rows = self._two_hump_rows()
        bdg_path = self._write_bedgraph(tmp, rows)
        bw_path = self._write_bigwig(tmp, rows)
        (tmp / "peaks.bed").write_text("chr1\t100\t2900\n")

        H.process(self._args(tmp, bdg_path).__class__(**{
            **vars(self._args(tmp, bdg_path)),
            "out": str(tmp / "bdg_out")
        }))
        H.process(self._args(tmp, bw_path).__class__(**{
            **vars(self._args(tmp, bw_path)),
            "out": str(tmp / "bw_out")
        }))

        bdg_nfr = (tmp / "bdg_out.nfr.bed").read_text()
        bw_nfr  = (tmp / "bw_out.nfr.bed").read_text()
        assert bdg_nfr == bw_nfr

        bdg_nuc = (tmp / "bdg_out.nuc.bed").read_text()
        bw_nuc  = (tmp / "bw_out.nuc.bed").read_text()
        assert bdg_nuc == bw_nuc

