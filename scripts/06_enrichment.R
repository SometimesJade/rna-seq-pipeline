# Step 6: GO over-representation, GO GSEA, Hallmark GSEA, PPP gene sets
source(here::here("scripts", "config.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggplot2)
})
set.seed(seed)

res_df <- readRDS(file.path(rds_dir, "res_df.rds"))
res    <- readRDS(file.path(rds_dir, "res.rds"))

save_plot <- function(p, name, w = 8, h = 7) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

#Over-representation analysis (ORA)
tested   <- res_df[!is.na(res_df$padj), ]
universe <- unique(tested$ensembl)

run_ora <- function(genes, label) {
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = padj_cutoff, qvalueCutoff = 0.05, readable = TRUE)
  write.csv(as.data.frame(ego),
            file.path(table_dir, paste0("GO_ORA_", label, ".csv")), row.names = FALSE)

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No significant GO terms for ", label, " genes; skipping plot")
    return(invisible(NULL))
  }
  ego_s <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  save_plot(dotplot(ego_s, showCategory = 15) +
              ggtitle(paste0("GO BP (ORA): ", label, " genes")),
            paste0("06_GO_ORA_", label))
  invisible(ego)
}

up_genes   <- unique(tested$ensembl[tested$diffexpressed == "Upregulated"])
down_genes <- unique(tested$ensembl[tested$diffexpressed == "Downregulated"])
ego_up   <- run_ora(up_genes,   "upregulated")
ego_down <- run_ora(down_genes, "downregulated")

#Ranked gene lists for GSEA (Wald statistic, unshrunken)
ranks_df <- data.frame(ensembl = sub("\\..*$", "", rownames(res)), stat = res$stat)
ranks_df <- ranks_df[!is.na(ranks_df$stat), ]
ranks_df <- ranks_df[order(ranks_df$stat, decreasing = TRUE), ]
ranks_df <- ranks_df[!duplicated(ranks_df$ensembl), ]
ranks    <- setNames(ranks_df$stat, ranks_df$ensembl)

sym <- mapIds(org.Hs.eg.db, names(ranks), "SYMBOL", "ENSEMBL", multiVals = "first")
ranks_sym <- ranks[!is.na(sym)]
names(ranks_sym) <- sym[!is.na(sym)]
ranks_sym <- ranks_sym[!duplicated(names(ranks_sym))]

plot_gsea <- function(gsea, title, name) {
  p <- dotplot(gsea, showCategory = 10, split = ".sign", label_format = 45) +
    facet_grid(. ~ .sign) + ggtitle(title)
  save_plot(p, name, w = 11, h = 11)
}

#GO GSEA
gsea_go <- gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                 ont = "BP", pvalueCutoff = padj_cutoff, seed = TRUE)
gsea_go <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
write.csv(as.data.frame(gsea_go), file.path(table_dir, "GSEA_GO_BP.csv"), row.names = FALSE)
plot_gsea(gsea_go, "GSEA: GO Biological Process", "07_GSEA_GO_BP")

#Hallmark GSEA
gmt_file <- list.files(here::here("metadata"), pattern = "^h\\.all.*symbols\\.gmt$",
                       full.names = TRUE)
if (length(gmt_file) == 1) {
  hallmark <- read.gmt(gmt_file)
  gsea_h <- GSEA(ranks_sym, TERM2GENE = hallmark,
                 pvalueCutoff = padj_cutoff, seed = TRUE)
  write.csv(as.data.frame(gsea_h), file.path(table_dir, "GSEA_Hallmark.csv"), row.names = FALSE)
  plot_gsea(gsea_h, "GSEA: MSigDB Hallmark", "08_GSEA_Hallmark")

  for (set in c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")) {
    if (set %in% gsea_h$ID) {
      png(file.path(fig_dir, paste0("09_", set, ".png")), width = 8, height = 6,
          units = "in", res = 300)
      print(gseaplot2(gsea_h, geneSetID = set, title = set))
      dev.off()
    }
  }
} else {
  message("Hallmark GMT not found in metadata/; skipping Hallmark GSEA")
}

#Pentose phosphate pathway branches
ppp_sets <- data.frame(
  term = c(rep("PPP_oxidative", 4), rep("PPP_nonoxidative", 5)),
  gene = c("G6PD", "PGLS", "PGD", "H6PD",
           "RPIA", "RPE", "TKT", "TKTL1", "TALDO1")
)
gsea_ppp <- GSEA(ranks_sym, TERM2GENE = ppp_sets, minGSSize = 3,
                 pvalueCutoff = 1, seed = TRUE)
write.csv(as.data.frame(gsea_ppp), file.path(table_dir, "GSEA_PPP_branches.csv"),
          row.names = FALSE)