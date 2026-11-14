# Base-R validation of processed input files; no model fitting.
local({
  read_input <- function(name, ...) read.csv(file.path("input", name), check.names = FALSE, ...)
  required <- c("amplicon_16S_asv_TZ_AU.csv", "amplicon_16S_tax_TZ_AU.csv",
                "meta_data_TZ_AU.csv", "metaproteome_TZ_AU.csv")
  stopifnot(all(file.exists(file.path("input", required))))
  meta <- read_input("meta_data_TZ_AU.csv")
  expected <- c("Mweka1", "Mweka2", "Mweka3", "Marangu", "Maua", "Umbwe", "Styria")
  stopifnot(nrow(meta) == 21L, !anyDuplicated(meta$id),
            setequal(unique(meta$group), expected), all(table(meta$group) == 3L))
  counts <- read_input("amplicon_16S_asv_TZ_AU.csv", row.names = 1)
  tax <- read_input("amplicon_16S_tax_TZ_AU.csv", row.names = 1)
  stopifnot(!anyDuplicated(rownames(counts)), setequal(rownames(counts), rownames(tax)),
            all(meta$id %in% colnames(counts)))
  x <- as.matrix(counts[, meta$id, drop = FALSE])
  stopifnot(is.numeric(x), all(is.finite(x)), all(x >= 0), all(x == round(x)),
            all(colSums(x) > 0))
  prot <- read_input("metaproteome_TZ_AU.csv")
  stopifnot(all(c("og", "subRole", "mainRole", "superKingdom", "Phylum", "Family", meta$id) %in% names(prot)))
  y <- as.matrix(prot[, meta$id, drop = FALSE])
  stopifnot(is.numeric(y), all(is.finite(y)), all(y >= 0), all(abs(colSums(y) - 1) < 1e-5))
  cat("Input checks passed: 21 thalli, 7 populations, 3 replicates each.\n")
})
