# Yungas Viral Community Statistics

R workflow used to perform ecological analyses of viral communities recovered from microbial mats in the Yungas rainforest.

This repository contains the code used to generate diversity statistics and figures for the study:

**"Uncovering phage diversity and functional potential in the Yungas rainforest through metagenomics."**

---

## Analyses included

The script performs:

• Counts per million (CPM) normalization  
• Alpha diversity metrics (Richness, Shannon, Simpson)  
• Rarefaction curves (iNEXT)  
• Rank–abundance distributions  
• Bray–Curtis dissimilarity  
• Principal Coordinates Analysis (PCoA)  
• Heatmap of the most abundant viral contigs  

---

## Input data

The script expects a tab-separated table:

contig_counts_3samples.tsv

with the structure:

sample contig_id contig_len mapped_reads

Example:

HMNa NODE_1234_length_5000_cov_4.2 5000 124
HMAb NODE_5678_length_3000_cov_3.1 3000 55
HMR NODE_9012_length_6200_cov_5.0 6200 210

---

## Generation of the input table

The table is generated prior to the R workflow by mapping metagenomic reads to viral contigs.

Typical workflow:

bowtie2-build viral_contigs.fasta index
bowtie2 -x index -1 reads_R1.fq.gz -2 reads_R2.fq.gz | samtools sort -o sample.bam
samtools idxstats sample.bam > sample.idxstats

The idxstats outputs are merged into a single table.

---

## Figures produced

The script generates the following figures:

Fig1_AlphaDiversity
Fig2_Rarefaction_iNEXT
Fig3_RankAbundance
Fig4_PCoA_BrayCurtis
FigS1_Heatmap_topContigs_logCPM

---

## Requirements

R packages:

tidyverse
vegan
ggrepel
iNEXT
ggplot2

---

## Running the analysis

Run in R:

source("viral_community_statistics.R")

---

## Author

**Hayde Saracho**  
PhD candidate in Biological Sciences  
Bioinformatics and microbial genomics
