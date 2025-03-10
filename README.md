# HisTrader

## Introduction

HisTrader is a **Perl-based tool** designed to identify **Nucleosome-Free Regions (NFRs)** from **ChIP-Seq data targeting histone modifications**. It takes a **BedGraph file** and a **peak file** (called from the BedGraph) as input and outputs **BED files** identifying both **NFRs** and **nucleosome-occupied regions (NORs)** within the defined peaks.

NFRs are detected using a combination of **moving averages** and **second-order differencing** to identify valley-like patterns within peaks.

- Running the software **without parameters** applies the **least stringent** settings, maximizing the number of predicted NFRs.
- **Stricter analysis** can be customized **case by case** using the provided parameters.

For more details, please refer to our **publication**: [https://doi.org/10.1101/2020.03.12.989228].

---

## Installation

Ensure you have **Perl** installed on your system. HisTrader requires standard Perl modules and runs directly from the command line.

---

## Usage

```bash
perl HISTRADER.pl --bedGraph ChIP.bedGraph --peaks ChIP.bed
```

### Required Parameters:

| Flag         | Description |
|-------------|-------------|
| `--bedGraph` | ChIP-Seq signal file in **BedGraph** format |
| `--peaks`    | **Broad peak** file in **BED** format |

### Optional Parameters:

| Flag          | Description |
|--------------|-------------|
| `--genome`   | Genome **FASTA** file (used to extract DNA sequences from NFRs) |
| `--trim`     | Trim extracted **FASTA** sequences (requires `--genome` and `--trimSize`) |
| `--trimSize` | Length (bp) of extracted sequences, centered on each NFR (requires `--trim` and `--genome`) |
| `--out`      | Output file prefix (default: **Histrader**) |
| `--method`   | Method for detecting NFRs: **MA** (Moving Average) / **DIFF** (Differencing) / **BOTH** (default) |
| `--step`     | Fixed **step size** (bp) for converting ChIP-Seq profiles (default: **25 bp**) |
| `--minSize`  | Minimum peak size for NFR calling (default: **500 bp**) |
| `--nucSize`  | Estimated **nucleosome size** (bp), used for moving average smoothing (default: **150 bp**) |
| `--mergeMulti` | Step multiplier for merging (default: **3×step** = **75 bp**) |
| `--maMulti`  | Moving average multiplier for slow-moving average (default: **3×nucleosome size**) |
| `--pMax`     | **Fraction of peak height** to use as a threshold for filtering (default: **0** = 0%) |
| `--filter`   | Filter NFRs larger than a specified size (default: **1000 bp**) |
| `--maxValley` | Print only NFR at the max peak region. If the NFR is filtered there will be no NFR for that peak (default: FALSE) |
| `--outBG`    | Output the fixed-step **ChIP-Seq signal** within broad peaks (in BedGraph format) |
| `--help`     | Display usage information |

---

## Example Command

```bash
perl HISTRADER.pl --bedGraph ChIP.bedGraph --peaks ChIP.bed --method BOTH --step 25 --nucSize 150 --pMax 0.1 --out output_prefix
```

---

## Citation

If you use **HisTrader** in your research, please cite:

> **[HisTrader: A Tool to Identify Nucleosome Free Regions from ChIP-Seq of Histone Post-Translational Modifications]**  
> [Efythios Kirbizakis, Yifei Yan, Ansley Gnanapragasam, Xiaoyang Zhang, and and Swneke Bailey]  
> [biorxiv, 2020]  
> [https://doi.org/10.1101/2020.03.12.989228]  

---

