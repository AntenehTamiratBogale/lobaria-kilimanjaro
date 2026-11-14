# Install dependencies in the active R library; run from the repository root.
cran <- c("vegan", "ggplot2", "ggrepel", "gtable", "statmod")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
bioc <- c("limma", "DESeq2")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
cat("Dependencies available.\n")
print(sessionInfo())
