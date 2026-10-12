# Step 4: import Salmon quantifications and run differential expression
source(here::here("scripts", "config.R"))
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ---- Sample metadata ----
samples <- read.csv(sample_sheet)
samples$condition <- factor(samples$condition, levels = c(ref_level, treat_level))
rownames(samples) <- samples$sample

files <- file.path(quant_dir, paste0(samples$sample, "_quant"), "quant.sf")
names(files) <- samples$sample
stopifnot("Missing quant.sf files" = all(file.exists(files)))

#Transcript-to-gene mapping from GENCODE headers
tx_names <- read.delim(files[1])$Name
parts <- strsplit(tx_names, "|", fixed = TRUE)
tx2gene <- data.frame(
  TXNAME = tx_names,
  GENEID = vapply(parts, `[`, character(1), 2),
  TXTYPE = vapply(parts, `[`, character(1), 8)
)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("TXNAME", "GENEID")])
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)

# Remove structural/small RNAs that poly-A libraries don't capture reliably
structural_types <- c("rRNA", "rRNA_pseudogene", "Mt_rRNA", "Mt_tRNA", "misc_RNA",
                      "snRNA", "snoRNA", "scaRNA", "ribozyme", "vault_RNA",
                      "scRNA", "sRNA")
structural_genes <- unique(tx2gene$GENEID[tx2gene$TXTYPE %in% structural_types])
dds <- dds[!rownames(dds) %in% structural_genes, ]
message("Removed ", sum(rownames(txi$counts) %in% structural_genes), " structural RNA genes")

#Import and build DESeq2 object
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("TXNAME", "GENEID")])
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)

structural_genes <- unique(tx2gene$GENEID[tx2gene$TXTYPE %in% structural_types])
n_struct <- sum(rownames(dds) %in% structural_genes)
dds <- dds[!rownames(dds) %in% structural_genes, ]
message("Removed ", n_struct, " structural RNA genes")


# Pre-filter low-count genes
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]
message("Genes kept after filtering: ", sum(keep), " of ", length(keep))

#Differential expression
dds <- DESeq(dds)
coef_name <- paste0("condition_", treat_level, "_vs_", ref_level)
stopifnot(coef_name %in% resultsNames(dds))

res    <- results(dds, name = coef_name, alpha = padj_cutoff)   # unshrunken, has 'stat' for GSEA
resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)

#Annotate and classify
res_df <- as.data.frame(resLFC)
res_df$gene_id <- rownames(res_df)
res_df$ensembl <- sub("\\..*$", "", res_df$gene_id)
res_df$symbol  <- mapIds(org.Hs.eg.db, res_df$ensembl, "SYMBOL", "ENSEMBL",
                         multiVals = "first")

sig <- !is.na(res_df$padj) & res_df$padj < padj_cutoff
res_df$diffexpressed <- "Not Significant"
res_df$diffexpressed[sig & res_df$log2FoldChange >  lfc_cutoff] <- "Upregulated"
res_df$diffexpressed[sig & res_df$log2FoldChange < -lfc_cutoff] <- "Downregulated"

res_df <- res_df[order(res_df$padj),
                 c("gene_id", "ensembl", "symbol", "baseMean",
                   "log2FoldChange", "lfcSE", "pvalue", "padj", "diffexpressed")]

#Variance-stabilised data for plots
vsd <- vst(dds, blind = FALSE)

#Save outputs
saveRDS(dds,    file.path(rds_dir, "dds.rds"))
saveRDS(vsd,    file.path(rds_dir, "vsd.rds"))
saveRDS(res,    file.path(rds_dir, "res.rds"))
saveRDS(res_df, file.path(rds_dir, "res_df.rds"))

write.csv(res_df, file.path(table_dir, "DE_amlodipine_vs_control.csv"), row.names = FALSE)
write.csv(counts(dds, normalized = TRUE),
          file.path(table_dir, "normalized_counts.csv"))

message("Upregulated: ",   sum(res_df$diffexpressed == "Upregulated"),
        " | Downregulated: ", sum(res_df$diffexpressed == "Downregulated"))