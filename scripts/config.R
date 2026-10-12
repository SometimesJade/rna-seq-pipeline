# Shared settings for all analysis steps
library(here)

# Paths
quant_dir    <- here("quants")
sample_sheet <- here("metadata", "samples.csv")
results_dir  <- here("results")
table_dir    <- file.path(results_dir, "tables")
fig_dir      <- file.path(results_dir, "figures")
rds_dir      <- file.path(results_dir, "rds")
for (d in c(table_dir, fig_dir, rds_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Design
ref_level   <- "Control"
treat_level <- "Amlodipine"

# Thresholds
padj_cutoff <- 0.05
lfc_cutoff  <- 1
min_count   <- 10
min_samples <- 3

# Plot colours and reproducibility
cond_colors <- c(Control = "grey40", Amlodipine = "#D55E00")
seed <- 42