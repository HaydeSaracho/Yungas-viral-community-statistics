# Yungas viral community statistics

This repository contains the scripts and processed data used to analyze viral community structure in microbial mats from hot springs of the Yungas rainforest (Jujuy - Argentina).

## Overview

Three microbial mat types were analyzed:

- HMAb: green mat
- HMNa: black mat
- HMR: pink mat

Each mat type was sampled in duplicate.

## Workflow

### 1. Mapping and abundance estimation
Reads from each replicate were mapped independently against the viral contigs corresponding to each mat type using Bowtie2 v2.5.2. Per-base coverage was calculated with SAMtools v1.19.2 depth including zero-coverage positions. Mean depth per contig was used as a proxy for viral population abundance.

Script:
- `scripts/01_map_and_compute_mean_depth.sh`

### 2. Ecological analysis
Processed mean-depth tables were used in R to calculate:

- alpha diversity (richness, Shannon, Simpson)
- Bray–Curtis dissimilarity
- PCoA
- exploratory PERMANOVA

Script:
- `scripts/02_ecological_analysis.R`

## Repository contents

- `data/coverage/`: processed abundance tables (`mean_depth.txt`)
- `data/metadata/`: sample metadata
- `scripts/`: reproducible analysis scripts
- `results/`: tables and figures generated from the R workflow

## Data availability

Raw sequencing reads are available at ENA under accession:

**PRJEB88976**

Intermediate and depth BAM files are also not included in this repository.

## Reproducibility

Downstream analyses can be reproduced from the processed coverage files by running:

```bash
Rscript scripts/02_ecological_analysis.R
```

## Full pipeline execution

To reproduce the full workflow from raw data:

```bash
bash scripts/01_map_and_compute_mean_depth.sh
Rscript scripts/02_ecological_analysis.R
```

## Software environment

A Conda environment file is provided:

```bash
conda env create -f environment.yml
conda activate yungas_virome
```
