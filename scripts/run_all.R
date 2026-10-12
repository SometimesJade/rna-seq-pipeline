# Run the full R analysis and record package versions
steps <- c("04_deseq2.R", "05_qc_plots.R", "06_enrichment.R")
for (s in steps) {
  message("\n==== Running ", s, " ====")
  source(here::here("scripts", s), local = new.env())
}
writeLines(capture.output(sessionInfo()), here::here("results", "sessionInfo.txt"))
message("\nDone. Results are in results/")