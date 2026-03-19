library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggrepel)

# Ruta a archivos mean depth
files <- list.files("coverage", pattern = "\\.mean_depth\\.txt$", full.names = TRUE)
files

# Leer archivos
read_mean_depth <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("contig", "mean_depth")
  
  sample_name <- sub("\\.mean_depth\\.txt$", "", basename(file))
  df$sample <- sample_name
  df
}

depth_long <- bind_rows(lapply(files, read_mean_depth))

# Formato ancho
depth_wide <- depth_long %>%
  pivot_wider(names_from = sample, values_from = mean_depth, values_fill = 0)

depth_mat <- as.data.frame(depth_wide)
rownames(depth_mat) <- depth_mat$contig
depth_mat <- depth_mat[, -1]

# Ordenar columnas
depth_mat <- depth_mat[, c("HMAb_1", "HMAb_2", "HMNa_1", "HMNa_2", "HMR_1", "HMR_2")]

# Metadata
metadata <- data.frame(
  sample = c("HMAb_1", "HMAb_2", "HMNa_1", "HMNa_2", "HMR_1", "HMR_2"),
  mat = c("HMAb", "HMAb", "HMNa", "HMNa", "HMR", "HMR"),
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$sample

# Normalización por tamaño de biblioteca
rel_abund <- sweep(depth_mat, 2, colSums(depth_mat), FUN = "/")

# Vegan usa muestras en filas
rel_abund_t <- t(rel_abund)

# Alpha diversity
richness <- rowSums(rel_abund_t > 0)
shannon <- diversity(rel_abund_t, index = "shannon")
simpson <- diversity(rel_abund_t, index = "simpson")

alpha_div <- data.frame(
  sample = rownames(rel_abund_t),
  richness = richness,
  shannon = shannon,
  simpson = simpson,
  mat = metadata$mat
)

print(alpha_div)
write.table(alpha_div, "alpha_diversity_mean_depth.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Bray-Curtis
bray <- vegdist(rel_abund_t, method = "bray")

# PCoA
pcoa <- cmdscale(bray, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  sample = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2],
  mat = metadata$mat
)

print(pcoa_df)
write.table(pcoa_df, "pcoa_bray_mean_depth.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Gráfico PCoA
p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = mat)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample), size = 3) +
  scale_color_manual(values = c(
    "HMAb" = "#1b5e20",   # verde oscuro
    "HMNa" = "#424242",   # gris oscuro (mejor que negro puro)
    "HMR"  = "#d81b60"    # rosa
  )) +
  theme_minimal() +
  labs(
    title = "PCoA based on Bray–Curtis dissimilarity",
    x = "PC1",
    y = "PC2",
    color = "Mat type"
  )

print(p)
ggsave("PCoA_bray_mean_depth.png", p, width = 7, height = 5, dpi = 300)

# PERMANOVA exploratoria
adonis_res <- adonis2(bray ~ mat, data = metadata)
print(adonis_res)
capture.output(adonis_res, file = "PERMANOVA_mean_depth.txt")