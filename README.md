# Yungas viral community statistics

This repository contains the scripts and processed data used to analyze viral community structure in microbial mats from hot springs of the Yungas rainforest (Jujuy, Argentina).

## Biological context

Three microbial mat types were sampled from hot springs:

| Code  | Mat color | Replicates |
|-------|-----------|------------|
| HMAb  | Green     | HMAb_1, HMAb_2 |
| HMNa  | Black     | HMNa_1, HMNa_2 |
| HMR   | Pink      | HMR_1, HMR_2   |

Each mat type represents **one metagenome sequenced twice** (two technical/biological replicates). Replicates were processed and mapped independently.

---

## Workflow

### 1. Mapping and abundance estimation

Reads from each replicate were mapped independently against the viral contigs corresponding to each mat type using **Bowtie2 v2.5.2**. Only mapped reads were retained (`-F 4`). Per-base coverage was calculated with **SAMtools v1.19.2** `depth`, including zero-coverage positions. Mean depth per contig was used as a proxy for viral population abundance.

Script: `scripts/01_map_and_compute_mean_depth.sh`

### 2. Ecological analysis

Processed mean-depth tables were used in R to calculate:

- Alpha diversity: richness, Shannon index, Simpson index
- Bray–Curtis dissimilarity matrix
- PCoA (Principal Coordinates Analysis)
- Exploratory PERMANOVA (`adonis2`)

Script: `scripts/02_ecological_analysis.R`

---

## Repository structure

```
.
├── data/
│   ├── coverage/          # Processed abundance tables (mean_depth.txt, one per replicate)
│   └── metadata/          # Sample metadata (samples_metadata.tsv)
├── scripts/
│   ├── 01_map_and_compute_mean_depth.sh
│   └── 02_ecological_analysis.R
├── results/
│   ├── tables/            # Alpha diversity, PCoA coordinates, PERMANOVA
│   └── figures/           # PCoA plot (PNG)
├── environment.yml        # Conda environment
└── .gitignore
```

---

## Key results

### Alpha diversity

| Sample  | Richness | Shannon | Simpson |
|---------|----------|---------|---------|
| HMAb_1  | 22       | 2.96    | 0.941   |
| HMAb_2  | 22       | 2.93    | 0.936   |
| HMNa_1  | 62       | 3.47    | 0.932   |
| HMNa_2  | 62       | 3.43    | 0.921   |
| HMR_1   | 8        | 1.94    | 0.842   |
| HMR_2   | 8        | 2.06    | 0.870   |

The black mat (HMNa) harbors the most diverse viral community (~62 viral populations), followed by the green mat (HMAb, ~22) and the pink mat (HMR, ~8). Replicates are highly consistent within each mat type.

### Beta diversity (PCoA + PERMANOVA)

PCoA based on Bray–Curtis dissimilarity shows clear separation between all three mat types, with replicates clustering tightly together. PERMANOVA (R² = 0.977) indicates that mat type explains ~98% of viral community variation.

> **Note on PERMANOVA power:** With only n=2 replicates per group, the minimum achievable p-value is ~0.083, making p<0.05 mathematically unattainable regardless of effect size. Results should therefore be interpreted primarily by effect size (R²) rather than by the p-value threshold.

---

## Reproducibility

### Requirements

Set up the Conda environment before running any scripts:

```bash
conda env create -f environment.yml
conda activate yungas_virome
```

### Run the full pipeline (from raw reads)

```bash
bash scripts/01_map_and_compute_mean_depth.sh
Rscript scripts/02_ecological_analysis.R
```

### Run only the ecological analysis (from processed coverage files)

```bash
Rscript scripts/02_ecological_analysis.R
```

Both scripts must be run **from the project root directory**.

---

## Software environment

| Tool       | Version  | Purpose                        |
|------------|----------|--------------------------------|
| Bowtie2    | 2.5.2    | Read mapping                   |
| SAMtools   | 1.19.2   | BAM processing and depth calc  |
| R          | ≥ 4.0    | Ecological analysis            |
| vegan      | ≥ 2.6    | Diversity metrics, PERMANOVA   |
| ggplot2    | ≥ 3.4    | Visualization                  |
| ggrepel    | ≥ 0.9    | Label repelling in plots       |
| dplyr/tidyr| ≥ 1.0    | Data manipulation              |

---

## Data availability

Raw sequencing reads are available at ENA under accession **PRJEB88976**.

Intermediate files (BAM, per-base depth) are not included in this repository (see `.gitignore`).
