library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggrepel)

# Reproducibility
set.seed(42)

# ----------------------------- Paths ---------------------------------------- #

# Paths relative to project root (run as: Rscript scripts/02_ecological_analysis.R)
cov_dir     <- "data/coverage"
out_tables  <- "results/tables"
out_figures <- "results/figures"

dir.create(out_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_figures, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- Load data ------------------------------------ #

files <- list.files(cov_dir, pattern = "\\.mean_depth\\.txt$", full.names = TRUE)

if (length(files) == 0) {
  stop("No mean_depth.txt files found in: ", cov_dir)
}

message("Found ", length(files), " coverage files:")
message(paste(" -", basename(files), collapse = "\n"))

read_mean_depth <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("contig", "mean_depth")
  df$sample <- sub("\\.mean_depth\\.txt$", "", basename(file))
  df
}

depth_long <- bind_rows(lapply(files, read_mean_depth))

# ----------------------------- Wide format ---------------------------------- #

depth_wide <- depth_long %>%
  pivot_wider(names_from = sample, values_from = mean_depth, values_fill = 0)

depth_mat <- as.data.frame(depth_wide)
rownames(depth_mat) <- depth_mat$contig
depth_mat <- depth_mat[, -1]

# Enforce consistent sample order
sample_order <- c("HMAb_1", "HMAb_2", "HMNa_1", "HMNa_2", "HMR_1", "HMR_2")
depth_mat <- depth_mat[, sample_order]

# ----------------------------- Metadata ------------------------------------- #

# NOTE: Each mat type represents one metagenome sampled twice (two replicates).
# Replicates are treated as independent observations in ordination and diversity
# calculations, but PERMANOVA power is limited by n=2 per group (see below).
metadata <- data.frame(
  sample = sample_order,
  mat    = c("HMAb", "HMAb", "HMNa", "HMNa", "HMR", "HMR"),
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$sample

# ----------------------------- Normalization -------------------------------- #

# Relative abundance (proportion of total library depth per sample)
rel_abund   <- sweep(depth_mat, 2, colSums(depth_mat), FUN = "/")
rel_abund_t <- t(rel_abund)   # vegan expects samples in rows

# ----------------------------- Alpha diversity ------------------------------ #

richness <- rowSums(rel_abund_t > 0)
shannon  <- diversity(rel_abund_t, index = "shannon")
simpson  <- diversity(rel_abund_t, index = "simpson")

alpha_div <- data.frame(
  sample   = rownames(rel_abund_t),
  richness = richness,
  shannon  = shannon,
  simpson  = simpson,
  mat      = metadata[rownames(rel_abund_t), "mat"]
)

print(alpha_div)
write.table(alpha_div,
            file      = file.path(out_tables, "alpha_diversity_mean_depth.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# ----------------------------- Beta diversity ------------------------------- #

bray <- vegdist(rel_abund_t, method = "bray")

# PCoA
pcoa    <- cmdscale(bray, eig = TRUE, k = 2)

# Percent variance explained per axis
eig_vals  <- pcoa$eig
pct_var   <- round(100 * eig_vals / sum(eig_vals[eig_vals > 0]), 1)
pc1_label <- paste0("PC1 (", pct_var[1], "% variance explained)")
pc2_label <- paste0("PC2 (", pct_var[2], "% variance explained)")

pcoa_df <- data.frame(
  sample = rownames(pcoa$points),
  PC1    = pcoa$points[, 1],
  PC2    = pcoa$points[, 2],
  mat    = metadata[rownames(pcoa$points), "mat"]
)

print(pcoa_df)
write.table(pcoa_df,
            file      = file.path(out_tables, "pcoa_bray_mean_depth.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# PCoA plot
p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = mat)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample), size = 3) +
  scale_color_manual(values = c(
    "HMAb" = "#1b5e20",   # dark green
    "HMNa" = "#424242",   # dark grey
    "HMR"  = "#d81b60"    # pink/rose
  )) +
  theme_minimal() +
  labs(
    title = "PCoA based on Bray–Curtis dissimilarity",
    x     = pc1_label,
    y     = pc2_label,
    color = "Mat type"
  )

print(p)
ggsave(file.path(out_figures, "PCoA_bray_mean_depth.png"),
       p, width = 7, height = 5, dpi = 300)

# ----------------------------- PERMANOVA ------------------------------------ #

# NOTE ON POWER: With n=2 replicates per group, the minimum achievable p-value
# under 719 permutations is ~0.083. A p-value < 0.05 is mathematically
# unattainable regardless of effect size. Results should therefore be interpreted
# by effect size (R²) rather than significance threshold.
adonis_res <- adonis2(bray ~ mat, data = metadata, permutations = 719)
print(adonis_res)

# Save with power note appended
permanova_out <- file.path(out_tables, "PERMANOVA_mean_depth.txt")
capture.output(adonis_res, file = permanova_out)
cat(
  "\nNOTE: With n=2 replicates per group, the minimum achievable p-value is ~0.083.",
  "\nInterpret by effect size (R²=", round(adonis_res$R2[1], 4), ") rather than",
  "\nthe p-value threshold of 0.05.\n",
  file  = permanova_out,
  append = TRUE
)

message("\nAll outputs written to:")
message("  Tables:  ", out_tables)
message("  Figures: ", out_figures)
