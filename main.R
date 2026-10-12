library(tximport)
library(readr)

#Quants Directory Location
dir <- here::here("quants")

#Samples
samples <- c("SRR40892062", "SRR40892063", "SRR40892064",
             "SRR40892065", "SRR40892066", "SRR40892067")

#Create and name file paths
files <- file.path(dir, paste0(samples, "_quant"), "quant.sf")
names(files) <- samples

#Finding Ensembl Gene ID from GENCODE header for importing
quant_data <- read.delim(files[1])
raw_transcript_names <- quant_data$Name
split_names <- strsplit(as.character(raw_transcript_names), "\\|")

#Create the dataframe for importing
tx2gene <- data.frame(
  TXNAME = raw_transcript_names,
  GENEID = sapply(split_names, `[`, 2)
)
#Importing with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

#Create Metadata dataframe
colData <- data.frame(
  row.names = samples,
  condition = factor(c(rep("Amlodipine", 3), rep("Control", 3)))
)

#Construct the DESeqDataSet
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "Control")

#Running DESeq2
dds <- DESeq(dds)

#Apply LFC Shrinkage
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_Amlodipine_vs_Control", type="apeglm")

#Variance Stabilization and PCA
library(ggplot2)

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Generate PCA Plot
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA Plot: Amlodipine vs Control") +
  labs(color = "Treatment")

#Generate the Volcano Plot
library(ggrepel)
res_df <- as.data.frame(resLFC)
res_df$gene <- rownames(res_df)

res_df$diffexpressed <- "Not Significant"
res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Upregulated"
res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Downregulated"
top <- head(res_df[order(res_df$padj), ], 15)

ggplot(res_df, aes(log2FoldChange, -log10(padj), col = diffexpressed)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = top, aes(label = symbol), color = "black", size = 3) +
  scale_color_manual(values = c(Downregulated = "blue",
                                "Not Significant" = "grey",
                                Upregulated = "red")) +
  theme_minimal()

#Go Enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

#Clean LFC
resLFC_clean <- resLFC[!is.na(resLFC$padj), ]
universe_gene_ids <- unique(gsub("\\..*", "", rownames(resLFC_clean)))

#Up-regulated genes
up_genes <- resLFC_clean[resLFC_clean$padj < 0.05 & resLFC_clean$log2FoldChange > 1, ]
up_gene_ids <- gsub("\\..*", "", rownames(up_genes))

#Down-regulated genes
down_genes <- resLFC_clean[resLFC_clean$padj < 0.05 & resLFC_clean$log2FoldChange < -1, ]
down_gene_ids <- gsub("\\..*", "", rownames(down_genes))

#Up-regulated GO-Enrichment
ego_up <- enrichGO(gene          = up_gene_ids,
                   universe      = universe_gene_ids,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

dotplot(ego_up, showCategory=15, title="Upregulated Pathways: Amlodipine")
ego_up_simple <- simplify(ego_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
dotplot(ego_up_simple, showCategory = 15)

#Down-regulated GO-enrichment
ego_down <- enrichGO(gene          = down_gene_ids,
                     universe      = universe_gene_ids,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

dotplot(ego_down, showCategory=15, title="Downregulated Pathways: Amlodipine")

#Down-regulated debug
length(down_gene_ids); length(up_gene_ids)
sum(down_gene_ids %in% universe_gene_ids)   # should equal length(down_gene_ids)

ego_down_all <- enrichGO(gene = down_gene_ids, universe = universe_gene_ids,
                         OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                         pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
head(as.data.frame(ego_down_all), 15)

#GSEA
res <- results(dds, name = "condition_Amlodipine_vs_Control")
ranks <- res$stat
names(ranks) <- gsub("\\..*", "", rownames(res))
ranks <- ranks[!is.na(ranks) & !duplicated(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)

gsea_bp <- gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                 ont = "BP", pvalueCutoff = 0.05, seed = TRUE)
dotplot(gsea_bp, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)

#GSEA MSigbdr
library(msigdbr)
sym <- mapIds(org.Hs.eg.db, names(ranks), "SYMBOL", "ENSEMBL")
ranks_sym <- ranks[!is.na(sym)]
names(ranks_sym) <- sym[!is.na(sym)]
ranks_sym <- ranks_sym[!duplicated(names(ranks_sym))]

hallmark <- read.gmt("metadata/h.all.v2026.1.Hs.symbols.gmt")
gsea_h <- GSEA(ranks_sym, TERM2GENE = hallmark, seed = TRUE)
dotplot(gsea_h, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)

#Custom gene set
ppp_sets <- data.frame(
  term = c(rep("PPP_oxidative", 4), rep("PPP_nonoxidative", 5)),
  gene = c("G6PD", "PGLS", "PGD", "H6PD",
           "RPIA", "RPE", "TKT", "TKTL1", "TALDO1")
)
gsea_ppp <- GSEA(ranks_sym, TERM2GENE = ppp_sets, minGSSize = 3,
                 pvalueCutoff = 1, seed = TRUE)
as.data.frame(gsea_ppp)[, c("ID", "setSize", "NES", "pvalue")]

#GSEA Lung breaching genes
gsea_df <- as.data.frame(gsea_bp)
gsea_df[grepl("lung|respiratory|organ growth", gsea_df$Description),
        c("Description", "setSize", "NES", "p.adjust", "core_enrichment")]

#Pentose Phospate Pathway
library(AnnotationDbi)
res_df$ensembl <- gsub("\\..*", "", res_df$gene)
res_df$symbol  <- mapIds(org.Hs.eg.db, res_df$ensembl, "SYMBOL", "ENSEMBL")

ppp <- c("G6PD", "PGLS", "PGD", "RPIA", "RPE", "TKT", "TKTL1", "TALDO1")
res_df[res_df$symbol %in% ppp, c("symbol", "log2FoldChange", "padj")]

#p53 target
p53_targets <- c("CDKN1A", "MDM2", "GDF15", "SESN1", "TP53I3", "BAX", "FAS", "RRM2B")
res_df[res_df$symbol %in% p53_targets, c("symbol", "log2FoldChange", "padj")]

library(pheatmap)
ids <- res_df$gene[res_df$symbol %in% p53_targets]
mat <- assay(vsd)[ids, ]
rownames(mat) <- res_df[ids, "symbol"]
mat <- mat - rowMeans(mat)          # center each gene
pheatmap(mat, annotation_col = colData, cluster_cols = FALSE,
         main = "p53 target genes (VST, centered)")