#!/usr/bin/env Rscript
# make_synthetic_data.R -- generate a synthetic methylation benchmark on its own
#
#   Rscript scripts/make_synthetic_data.R [out_dir] [n_cpg] [n_signal]
#
# Defaults: data/synthetic/custom, 20000 CpGs, 100 planted signal CpGs with
# |delta-beta| 0.02-0.08, discovery n = 186, validation n = 90 (the design of
# the original study). See R/synthetic_data.R for the remaining knobs.

args <- commandArgs(trailingOnly = TRUE)
full <- commandArgs(trailingOnly = FALSE)
f    <- sub("^--file=", "", full[grepl("^--file=", full)])
root <- normalizePath(file.path(if (length(f)) dirname(f[1]) else getwd(), ".."))
if (!file.exists(file.path(root, "R", "utils.R"))) root <- getwd()
source(file.path(root, "R", "utils.R"))
source(file.path(root, "R", "synthetic_data.R"))

out_dir  <- if (length(args) >= 1) args[1] else file.path(root, "data", "synthetic", "custom")
n_cpg    <- if (length(args) >= 2) as.integer(args[2]) else 20000L
n_signal <- if (length(args) >= 3) as.integer(args[3]) else 100L

make_synthetic_methylation(out_dir = out_dir, n_cpg = n_cpg, n_signal = n_signal)
