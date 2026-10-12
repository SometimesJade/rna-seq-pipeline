# Step 5: QC and results figures
source(here::here("scripts", "config.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

vsd    <- readRDS(file.path(rds_dir, "vsd.rds"))
res_df <- readRDS(file.path(rds_dir, "res_df.rds"))

save_plot <- function(p, name, w = 7, h = 6) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
col_ann    <- data.frame(condition = vsd$condition, row.names = colnames(vsd))
ann_colors <- list(condition = cond_colors)

#PCA
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))

p_pca <- ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = name), size = 3, show.legend = FALSE) +
  scale_color_manual(values = cond_colors) +
  labs(x = paste0("PC1: ", pv[1], "% variance"),
       y = paste0("PC2: ", pv[2], "% variance"),
       color = "Treatment", title = "PCA: Amlodipine vs Control") +
  coord_fixed() +
  theme_bw()
save_plot(p_pca, "01_pca")

#Sample-to-sample distances
d <- dist(t(assay(vsd)))
pheatmap(as.matrix(d),
         clustering_distance_rows = d, clustering_distance_cols = d,
         annotation_col = col_ann, annotation_colors = ann_colors,
         main = "Sample-to-sample distances (VST)",
         filename = file.path(fig_dir, "02_sample_distances.png"),
         width = 7, height = 6)

#MA plot
de_colors <- c(Downregulated = "blue", "Not Significant" = "grey70", Upregulated = "red")

p_ma <- ggplot(res_df, aes(log10(baseMean), log2FoldChange, color = diffexpressed)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = de_colors) +
  labs(x = "log10 mean normalized count", y = "log2 fold change (apeglm)",
       color = NULL, title = "MA plot") +
  theme_minimal()
save_plot(p_ma, "03_ma_plot")

#Volcano plot with top genes labelled
plot_df <- res_df[!is.na(res_df$padj), ]
plot_df$label <- ifelse(is.na(plot_df$symbol), plot_df$ensembl, plot_df$symbol)
top <- head(plot_df[plot_df$diffexpressed != "Not Significant", ], 15)

p_volcano <- ggplot(plot_df, aes(log2FoldChange, -log10(padj), color = diffexpressed)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = top, aes(label = label), color = "black",
                  size = 3, max.overlaps = Inf) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  scale_color_manual(values = de_colors) +
  labs(x = "log2 fold change", y = "-log10(adjusted p-value)",
       color = NULL, title = "Volcano plot: Amlodipine vs Control") +
  theme_minimal()
save_plot(p_volcano, "04_volcano")

#Heatmap of genes of interest
genes_of_interest <- list(
  "p53 targets" = c("CDKN1A", "MDM2", "GDF15", "BAX", "FAS", "TP53I3", "SESN1", "RRM2B"),
  "Cell cycle"  = c("E2F1", "MKI67", "CCNB1", "PLK1"),
  "PPP"         = c("G6PD", "PGLS", "PGD", "RPIA", "RPE", "TKT", "TKTL1", "TALDO1")
)
goi <- stack(genes_of_interest)
goi <- goi[goi$values %in% res_df$symbol, ]
goi$ind <- droplevels(goi$ind)

ids <- res_df$gene_id[match(goi$values, res_df$symbol)]
mat <- assay(vsd)[ids, order(vsd$condition)]
mat <- mat - rowMeans(mat)
rownames(mat) <- goi$values

pheatmap(mat,
         cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_row = head(cumsum(table(goi$ind)), -1),
         annotation_col = col_ann, annotation_colors = ann_colors,
         annotation_row = data.frame(group = goi$ind, row.names = goi$values),
         main = "Selected genes (VST, centered)",
         filename = file.path(fig_dir, "05_genes_of_interest_heatmap.png"),
         width = 7, height = 8)