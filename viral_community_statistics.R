# ==============================================================================
# Viral community statistics – Yungas microbial mats
#
# Manuscript:
#   "Uncovering phage diversity and functional potential in Yungas rainforest
#    through metagenomics"
#
# Purpose:
#   Statistical analyses and figure generation for viral community ecology using
#   mapped-read counts per viral contig (treated as viral populations).
#
# Analyses included:
#   - CPM normalization (counts per million mapped reads to viral contigs)
#   - Alpha diversity: Observed richness, Shannon, Simpson
#   - Rarefaction (iNEXT; q=0; abundance-based)
#   - Rank–abundance curves (Whittaker plot)
#   - Bray–Curtis dissimilarity + PCoA
#   - Heatmap of top-N contigs (log10(CPM+1)) [supplementary]
#
# Input:
#   contig_counts_3samples.tsv with columns:
#     sample, contig_id, contig_len, mapped_reads
#
# Expected sample labels (either set works):
#   - Black / Green / Pink
#   - HMNa / HMAb / HMR
#
# Outputs:
#   - Table_alpha_diversity.tsv
#   - Fig1_AlphaDiversity.(png/pdf)
#   - Fig2_Rarefaction_iNEXT.(png/pdf)
#   - Fig3_RankAbundance.(png/pdf)
#   - Fig4_PCoA_BrayCurtis.(png/pdf)
#   - FigS1_Heatmap_topContigs_logCPM.(png/pdf)
#   - community_stats_objects.rds
#
# Notes:
#   - This workflow is descriptive (n=1 per mat). No hypothesis testing is done.
#   - Mapping/count generation is performed outside R (see below).
# ==============================================================================


# ==============================================================================
# Input file provenance
#
# The input table `contig_counts_3samples.tsv` is generated prior to this script
# by mapping paired-end metagenomic reads from each microbial mat to the set of
# viral contigs assembled for that same sample.
#
# Summary workflow (performed outside R):
#
# 1) Build Bowtie2 index from viral contigs for each mat/sample
#      bowtie2-build <viral_contigs.fasta> <index_prefix>
#
# 2) Map paired-end reads to viral contigs and generate sorted BAM
#      bowtie2 -x <index_prefix> -1 <reads_R1.fq.gz> -2 <reads_R2.fq.gz> \
#        --threads <N> | samtools view -bS - | samtools sort -o <sample>.bam
#
# 3) Index BAM and extract per-contig counts
#      samtools index <sample>.bam
#      samtools idxstats <sample>.bam > <sample>.idxstats.tsv
#
# 4) Merge per-sample idxstats outputs into a single TSV table with columns:
#      sample, contig_id, contig_len, mapped_reads
#
# Note: This R script performs statistical analyses and figure generation only.
# ==============================================================================


# ---------
# Packages
# ---------
pkgs <- c("tidyverse", "vegan", "ggrepel", "iNEXT", "ggplot2")
missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(missing) > 0) {
  stop(
    "Missing R packages: ", paste(missing, collapse = ", "),
    "\nInstall them with:\n  install.packages(c(",
    paste0('"', missing, '"', collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggrepel)
  library(iNEXT)
  library(ggplot2)
})


# ---------------------------
# User-configurable settings
# ---------------------------
infile  <- "contig_counts_3samples.tsv"
out_dir <- "."          # e.g., "outputs"
top_n   <- 30           # heatmap top-N contigs
dpi_png <- 600

# Consistent palette across all figures
mat_colors <- c(
  HMNa = "#2c3e50",  # black mat
  HMAb = "#27ae60",  # green mat
  HMR  = "#e84393"   # pink mat
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ------------------------
# Read input and validate
# ------------------------
if (!file.exists(infile)) stop("Input file not found: ", infile)

df <- readr::read_tsv(infile, show_col_types = FALSE)

required_cols <- c("sample", "contig_id", "contig_len", "mapped_reads")
if (!all(required_cols %in% colnames(df))) {
  stop("Input file must contain columns: ", paste(required_cols, collapse = ", "))
}

# Map sample labels to manuscript abbreviations
df <- df %>%
  mutate(sample = recode(sample,
                         "Black" = "HMNa",
                         "Green" = "HMAb",
                         "Pink"  = "HMR",
                         .default = sample)) %>%
  mutate(sample = factor(sample, levels = c("HMNa", "HMAb", "HMR")))

if (any(is.na(df$sample))) {
  stop("Some samples could not be mapped to HMNa/HMAb/HMR. Check the 'sample' column.")
}


# ---------------------
# Build matrices + CPM
# ---------------------
df_sum <- df %>%
  group_by(sample, contig_id) %>%
  summarise(mapped_reads = sum(mapped_reads), .groups = "drop")

mat <- df_sum %>%
  pivot_wider(names_from = contig_id,
              values_from = mapped_reads,
              values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()

mat <- mat[c("HMNa", "HMAb", "HMR"), , drop = FALSE]

# CPM normalization:
mat_cpm <- sweep(mat, 1, rowSums(mat), "/") * 1e6


# ----------------------
# Alpha diversity table
# ----------------------
alpha <- tibble(
  sample   = rownames(mat_cpm),
  richness = vegan::specnumber(mat_cpm),
  shannon  = vegan::diversity(mat_cpm, "shannon"),
  simpson  = vegan::diversity(mat_cpm, "simpson"),
  total_mapped_reads_to_viral_contigs = rowSums(mat)
) %>%
  mutate(sample = factor(sample, levels = c("HMNa", "HMAb", "HMR"))) %>%
  arrange(sample)

readr::write_tsv(alpha, file.path(out_dir, "Table_alpha_diversity.tsv"))


# --------------------------
# Figure 1: Alpha diversity
# --------------------------
alpha_long <- alpha %>%
  select(sample, richness, shannon, simpson) %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric,
                    richness = "Observed richness",
                    shannon  = "Shannon index",
                    simpson  = "Simpson index"),
    metric = factor(metric, levels = c("Observed richness", "Shannon index", "Simpson index"))
  )

p_alpha <- ggplot(alpha_long, aes(sample, value, fill = sample)) +
  geom_col(width = 0.7) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = mat_colors) +
  theme_classic(base_size = 14) +
  labs(title = "Alpha diversity across microbial mats", x = NULL, y = NULL) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "plain")
  )

ggsave(file.path(out_dir, "Fig1_AlphaDiversity.png"), p_alpha,
       width = 8, height = 3.2, dpi = dpi_png, bg = "white")
ggsave(file.path(out_dir, "Fig1_AlphaDiversity.pdf"), p_alpha,
       width = 8, height = 3.2, bg = "white")


# ----------------------------------
# Figure 2: Rarefaction (iNEXT; q=0)
# ----------------------------------
x_ab <- lapply(rownames(mat), function(s) as.numeric(mat[s, ]))
names(x_ab) <- rownames(mat)

rare_ab <- iNEXT::iNEXT(x_ab, q = 0, datatype = "abundance")

p_rare <- iNEXT::ggiNEXT(rare_ab, type = 1) +
  scale_color_manual(values = mat_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "Rarefaction curves of viral populations",
    x = "Sequencing depth (reads mapped to viral contigs)",
    y = "Estimated richness (q=0)"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Fig2_Rarefaction_iNEXT.png"), p_rare,
       width = 6.5, height = 5, dpi = dpi_png, bg = "white")
ggsave(file.path(out_dir, "Fig2_Rarefaction_iNEXT.pdf"), p_rare,
       width = 6.5, height = 5, bg = "white")


# -------------------------
# Figure 3: Rank–abundance
# -------------------------
rank_df <- as.data.frame(mat_cpm) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "contig_id", values_to = "cpm") %>%
  filter(cpm > 0) %>%
  group_by(sample) %>%
  arrange(desc(cpm), .by_group = TRUE) %>%
  mutate(rank = row_number(),
         rel  = cpm / sum(cpm)) %>%
  ungroup() %>%
  mutate(sample = factor(sample, levels = c("HMNa", "HMAb", "HMR")))

p_rank <- ggplot(rank_df, aes(rank, rel, color = sample)) +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  scale_color_manual(values = mat_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "Rank–abundance distribution of viral populations",
    x = "Rank (most abundant → least abundant)",
    y = "Relative abundance (log10)"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Fig3_RankAbundance.png"), p_rank,
       width = 6.5, height = 5, dpi = dpi_png, bg = "white")
ggsave(file.path(out_dir, "Fig3_RankAbundance.pdf"), p_rank,
       width = 6.5, height = 5, bg = "white")


# -----------------------------
# Figure 4: Bray–Curtis + PCoA
# -----------------------------
dist_bc <- vegan::vegdist(mat_cpm, method = "bray")

pcoa <- cmdscale(dist_bc, k = 2, eig = TRUE)
var_expl <- round(100 * pcoa$eig / sum(pcoa$eig), 1)

pcoa_df <- data.frame(
  sample = rownames(mat_cpm),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  mutate(sample = factor(sample, levels = c("HMNa", "HMAb", "HMR")))

p_pcoa <- ggplot(pcoa_df, aes(PC1, PC2, color = sample)) +
  geom_point(size = 5) +
  geom_text_repel(
    aes(label = sample),
    size = 6,
    box.padding = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(values = mat_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCoA of viral community composition",
    x = paste0("PC1 (", var_expl[1], "%)"),
    y = paste0("PC2 (", var_expl[2], "%)")
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    plot.margin = margin(10, 20, 10, 20)
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15))

ggsave(file.path(out_dir, "Fig4_PCoA_BrayCurtis.png"), p_pcoa,
       width = 6, height = 5, dpi = dpi_png, bg = "white")
ggsave(file.path(out_dir, "Fig4_PCoA_BrayCurtis.pdf"), p_pcoa,
       width = 6, height = 5, bg = "white")


# ---------------------------------------------------
# Supplementary Figure: Heatmap (top-N; log10(CPM+1))
# ---------------------------------------------------
mat_cpm_df <- as.data.frame(mat_cpm) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "contig_id", values_to = "cpm")

top_contigs <- mat_cpm_df %>%
  group_by(contig_id) %>%
  summarise(total_cpm = sum(cpm), .groups = "drop") %>%
  arrange(desc(total_cpm)) %>%
  slice_head(n = top_n) %>%
  pull(contig_id)

short_lab <- function(x) sub("^(NODE_[0-9]+).*", "\\1", x)

p_hm <- mat_cpm_df %>%
  filter(contig_id %in% top_contigs) %>%
  mutate(
    sample = factor(sample, levels = c("HMNa", "HMAb", "HMR")),
    contig_id = factor(contig_id, levels = top_contigs),
    contig_short = factor(short_lab(as.character(contig_id)),
                          levels = short_lab(top_contigs))
  ) %>%
  ggplot(aes(contig_short, sample, fill = log10(cpm + 1))) +
  geom_tile() +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  labs(
    title = paste0("Viral contig abundance heatmap (top ", top_n, ", log10(CPM+1))"),
    fill = "log10(CPM+1)"
  )

ggsave(file.path(out_dir, "FigS1_Heatmap_topContigs_logCPM.png"), p_hm,
       width = 12, height = 3.6, dpi = dpi_png, bg = "white")
ggsave(file.path(out_dir, "FigS1_Heatmap_topContigs_logCPM.pdf"), p_hm,
       width = 12, height = 3.6, bg = "white")


# -----------------------
# Save objects for reuse
# -----------------------
saveRDS(
  list(alpha = alpha, mat = mat, mat_cpm = mat_cpm, dist_bc = dist_bc, pcoa_df = pcoa_df),
  file = file.path(out_dir, "community_stats_objects.rds")
)


# -----------------------------------
# Console outputs (copy into Results)
# -----------------------------------
message("Done. Outputs written to: ", normalizePath(out_dir))
print(alpha)
print(as.matrix(dist_bc))

# Optional (recommended for reproducibility in your own runs)
# writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
