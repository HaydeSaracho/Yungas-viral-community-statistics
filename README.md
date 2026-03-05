# Yungas-viral-community-statistics
R script for viral community diversity analyses from Yungas microbial mats virome study

# Viral community statistics – Yungas microbial mats

This repository contains the R script used to perform statistical analyses and generate figures for the manuscript:

**Uncovering phage diversity and functional potential in Yungas rainforest through metagenomics**

The script reproduces the statistical analyses used to characterize viral community diversity across microbial mats from the Yungas rainforest.

---

## Analyses included

The workflow performs the following analyses:

- Alpha diversity (Observed richness, Shannon, Simpson)
- Rarefaction curves using **iNEXT**
- Rank–abundance distributions
- Bray–Curtis dissimilarity
- Principal Coordinates Analysis (PCoA)
- Heatmap of viral contig abundances

---

## Software requirements

The analysis was performed in **R** using the following packages:
tidyverse
vegan
iNEXT
ggrepel
ggplot2

---

## Input data format

The script expects a tab-separated table with the following columns:

| Column | Description |
|------|------|
| sample | microbial mat identifier |
| contig_id | viral contig identifier |
| contig_len | contig length |
| mapped_reads | number of reads mapped to the contig |

---

## Notes

This repository provides the statistical workflow used to generate the figures presented in the manuscript.

Raw sequencing data and assemblies are available through the corresponding study resources.
