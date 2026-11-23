#!/usr/bin/env Rscript
# Portable entry point for the study's original R implementation.
args <- commandArgs(trailingOnly = TRUE)
entry <- grep("^--file=", commandArgs(), value = TRUE)
if (length(entry) != 1L) stop("Run reproduce.R with Rscript.")
base <- dirname(normalizePath(sub("^--file=", "", entry)))
source(file.path(base, "raw_cusum_full.R"), local = TRUE, encoding = "UTF-8")
data_dir <- file.path(base, "data")
if (length(args) == 5L && args[1] == "chunk") {
    ids <- suppressWarnings(as.numeric(args[2:4]))
    if (anyNA(ids) || any(!is.finite(ids)) || any(ids != floor(ids))) {
        stop("Use integer cell, start and count. Count is at most 10.")
    }
    run_chunk("full", as.integer(ids[1]), as.integer(ids[2]), as.integer(ids[3]),
              file.path(data_dir, "seeds.rds"), args[5])
} else if (length(args) == 2L && args[1] == "tables") {
    if (file.exists(args[2])) stop("Choose a new output directory.")
    rows <- read.csv(file.path(data_dir, "trials.csv"), check.names = FALSE)
    stopifnot(nrow(rows) == 24000L, length(unique(rows$cell)) == 12L,
              all(table(rows$cell) == 2000L), !anyDuplicated(rows[c("cell", "replication")]))
    dir.create(args[2], recursive = TRUE)
    summarize_results(rows, get_design("full"), args[2])
} else {
    stop(paste("Usage: Rscript reproduce.R chunk CELL START COUNT NEW_DIR",
               "or Rscript reproduce.R tables NEW_DIR."))
}
