# HisTrader

## Introduction

HisTrader is a tool designed to identify **Nucleosome-Free Regions (NFRs)** from **ChIP-Seq data targeting histone modifications**. It takes a signal file (bedGraph or bigWig) and a broad peak file as input, and outputs BED files identifying **NFRs** and **nucleosome-occupied regions (NORs)** within those peaks.

NFRs are detected using a combination of **moving averages** and **second-order differencing** to identify valley-like patterns within peaks.

- Running without optional parameters applies the **least stringent** settings, maximising the number of predicted NFRs.
- Stricter analysis can be customised using the parameters below.

For more details, please refer to our publication: https://doi.org/10.1101/2020.03.12.989228

---

## Installation

Requires Python 3.8+ and the following packages:

```bash
pip install histrader

# or

pip install numpy pandas pyranges pyBigWig pytest
pip install .
```

---

## Usage

```bash
python histrader.py --signal ChIP.bedGraph/bigWig --peaks ChIP.bed
```

### Required Parameters

| Flag | Description |
|------|-------------|
| `--signal` | ChIP-Seq signal file in **bedGraph** or **bigWig** format |
| `--peaks` | Broad peak file in **BED** format |

### Optional Parameters

| Flag | Description | Default |
|------|-------------|---------|
| `--genome` | Genome FASTA file (used to extract DNA sequences from NFRs) | — |
| `--trim` | Trim extracted FASTA sequences (requires `--genome` and `--trimSize`) | off |
| `--trimSize` | Length (bp) of trimmed sequences, centred on each NFR | 100 |
| `--out` | Output file prefix | `Histrader` |
| `--method` | NFR detection method: `MA`, `DIFF`, or `BOTH` | `BOTH` |
| `--step` | Fixed step size (bp) for the converted signal | 25 |
| `--minSize` | Minimum peak size (bp) for NFR calling | 500 |
| `--nucSize` | Estimated nucleosome size (bp); should be divisible by `--step` | 150 |
| `--mergeMulti` | Step multiplier for merging (merge distance = mergeMulti × step) | 3 |
| `--maMulti` | Nucleosome multiplier for the slow moving average | 3 |
| `--pMax` | Fraction of peak max signal used as a zeroing threshold (0 = off) | 0.0 |
| `--filter` | Discard NFRs wider than this value (bp) | 1000 |
| `--maxValley` | Keep only the NFR at the max peak region per peak | off |
| `--useDifferential` | With `--maxValley`, use differential-based NFR selection | off |
| `--randValley` | Keep only one random NFR per peak | off |
| `--seed` | Random seed for reproducible `--randValley` output | — |
| `--outBG` | Also output the fixed-step signal within peaks (bedGraph) | off |

---

## Output Files

| File | Contents |
|------|----------|
| `<prefix>.nfr.bed` | Nucleosome-Free Regions |
| `<prefix>.nuc.bed` | Nucleosome-Occupied Regions |
| `<prefix>.missing.bed` | Peaks where no NFR was detected |
| `<prefix>.nfr.fa` | NFR sequences (only with `--genome`) |
| `<prefix>.nuc.fa` | Nucleosome sequences (only with `--genome`) |
| `<prefix>.bedGraph` | Fixed-step signal track (only with `--outBG`) |

---

## Example

```bash
python histrader.py \
  --signal  ChIP.bedGraph \
  --peaks   ChIP.bed \
  --method  BOTH \
  --step    25 \
  --nucSize 150 \
  --pMax    0.1 \
  --out     output_prefix
```

---

## Test Data

Test data for a single H3K27ac peak on chr11 is provided in `TEST_DATA/`:

```bash
python histrader.py \
  --signal TEST_DATA/test.region.histrader.chr11_12286115_12289133.H3K27AC.bdg \
  --peaks  TEST_DATA/test.region.histrader.chr11_12286115_12289133.H3K27AC.broadPeak \
  --out    test_output
```

## Tests

```bash
pytest test_histrader.py -v
```

---

## Citation

If you use HisTrader in your research, please cite:

> **HisTrader: A Tool to Identify Nucleosome Free Regions from ChIP-Seq of Histone Post-Translational Modifications**
> Eftyhios Kirbizakis, Yifei Yan, Ansley Gnanapragasam, Xiaoyang Zhang, and Swneke Bailey
> bioRxiv, 2020
> https://doi.org/10.1101/2020.03.12.989228
