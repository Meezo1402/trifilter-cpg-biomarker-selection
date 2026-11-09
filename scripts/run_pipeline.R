#!/usr/bin/env Rscript
# run_pipeline.R -- end-to-end driver for the permutational tri-filter pipeline
#
# Usage
#   Rscript scripts/run_pipeline.R --demo                        synthetic data, moderate settings
#   Rscript scripts/run_pipeline.R --smoke                       tiny settings for CI (~1 min)
#   Rscript scripts/run_pipeline.R --config config/my_study.R    your own cohorts
#   ... --out some/dir                                           override the output directory
#
# Stages
#   1. Tri-filter on the discovery cohort  -> candidate CpG pool
#   2. Permutation null from random panels -> what "chance" looks like on validation
#   3. RFE with RF / linear SVM / radial SVM -> selected panels
#   4. Validation scoring + empirical p-values against the null -> summary

`%||%`  <- function(a, b) if (is.null(a)) b else a
args     <- commandArgs(trailingOnly = TRUE)
has_flag <- function(flag) flag %in% args
get_arg  <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }

find_root <- function() {
  full  <- commandArgs(trailingOnly = FALSE)
  f     <- sub("^--file=", "", full[grepl("^--file=", full)])
  start <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
  for (cand in c(start, dirname(start))) {
    if (file.exists(file.path(cand, "R", "utils.R"))) return(normalizePath(cand))
  }
  stop("Could not locate the repository root; run from the repo root or via Rscript scripts/run_pipeline.R")
}
root <- find_root()
for (f in c("utils.R", "config.R", "01_trifilter.R", "02_permutation_classifier.R",
            "03_rfe_benchmark.R", "04_evaluate.R", "synthetic_data.R")) {
  source(file.path(root, "R", f))
}
check_packages(c("caret", "randomForest", "kernlab", "pROC", "e1071"))

# ---- configuration -----------------------------------------------------------
if (has_flag("--config")) {
  env <- new.env()
  sys.source(normalizePath(get_arg("--config")), envir = env)
  config <- env$config
  if (is.null(config)) stop("The config file must assign a list to an object named `config`")
  config$mode <- "custom"
} else {
  config <- default_config(root, if (has_flag("--smoke")) "smoke" else "demo")
}
out_dir <- get_arg("--out", config$out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(config$seed)
log_msg("mode: %s | output: %s", config$mode, out_dir)

# ---- data ----------------------------------------------------------------------
if (isTRUE(config$data$synthetic) && !file.exists(config$data$discovery_beta)) {
  log_msg("stage 0: generating synthetic cohorts")
  do.call(make_synthetic_methylation,
          c(list(out_dir = dirname(config$data$discovery_beta),
                 positive = config$data$positive, negative = config$data$negative,
                 seed = config$seed), config$synthetic))
}
log_msg("stage 0: loading cohorts")
disc <- with(config$data, load_cohort(discovery_beta,  discovery_labels,  sample_col, outcome_col, positive, negative))
vali <- with(config$data, load_cohort(validation_beta, validation_labels, sample_col, outcome_col, positive, negative))
shared <- intersect(colnames(disc$x), colnames(vali$x))
disc$x <- disc$x[, shared, drop = FALSE]
vali$x <- vali$x[, shared, drop = FALSE]
log_msg("  discovery n=%d (%s), validation n=%d (%s), %d CpGs shared",
        nrow(disc$x), paste(sprintf("%s=%d", levels(disc$y), table(disc$y)), collapse = ", "),
        nrow(vali$x), paste(sprintf("%s=%d", levels(vali$y), table(vali$y)), collapse = ", "),
        length(shared))

# ---- stage 1: tri-filter (discovery only) ----------------------------------------
log_msg("stage 1: tri-filter feature pre-selection")
tf <- do.call(trifilter_select, c(list(x = disc$x, y = disc$y), config$trifilter))
write.csv(tf$table, file.path(out_dir, "trifilter_table.csv"), row.names = FALSE)
writeLines(tf$features, file.path(out_dir, "trifilter_selected_cpgs.txt"))
pool <- tf$features
if (length(pool) < config$permutation$panel_size) {
  stop("Only ", length(pool), " CpGs survived the tri-filter; relax the thresholds or lower panel_size")
}
x_train <- disc$x[, pool, drop = FALSE]
x_valid <- vali$x[, pool, drop = FALSE]

# ---- stage 2: permutation null -----------------------------------------------------
log_msg("stage 2: permutation null distribution")
null <- do.call(permutation_null,
                c(list(x_train = x_train, y_train = disc$y, x_valid = x_valid, y_valid = vali$y,
                       feature_pool = pool, seed = config$seed, n_cores = config$n_cores),
                  config$permutation))
write.csv(null, file.path(out_dir, "permutation_null.csv"), row.names = FALSE)

# ---- stage 3: RFE benchmark -----------------------------------------------------------
log_msg("stage 3: recursive feature elimination")
rfe_res <- do.call(rfe_benchmark,
                   c(list(x_train = x_train, y_train = disc$y, x_valid = x_valid, y_valid = vali$y,
                          seed = config$seed), config$rfe))
for (nm in names(rfe_res)) {
  r <- rfe_res[[nm]]
  write.csv(r$profile, file.path(out_dir, sprintf("rfe_%s_profile.csv", nm)), row.names = FALSE)
  writeLines(r$selected, file.path(out_dir, sprintf("rfe_%s_selected_cpgs.txt", nm)))
  write.csv(r$validation$predictions, file.path(out_dir, sprintf("rfe_%s_validation_predictions.csv", nm)), row.names = FALSE)
}

# ---- stage 4: summary ------------------------------------------------------------------
log_msg("stage 4: summary")
truth <- if (!is.null(config$data$truth_signal) && file.exists(config$data$truth_signal)) readLines(config$data$truth_signal) else NULL

summary_tbl <- do.call(rbind, lapply(rfe_res, function(r) {
  v <- r$validation
  data.frame(learner = r$learner, n_selected = length(r$selected),
             cv_roc = r$cv_roc, valid_auc = v$auc, valid_accuracy = v$accuracy,
             valid_balanced_accuracy = v$balanced_accuracy,
             valid_sensitivity = v$sensitivity, valid_specificity = v$specificity,
             empirical_p_auc = empirical_p(null$valid_auc, v$auc),
             empirical_p_accuracy = empirical_p(null$valid_accuracy, v$accuracy),
             planted_signal_recovered = if (is.null(truth)) NA_integer_ else sum(r$selected %in% truth),
             stringsAsFactors = FALSE)
}))
write.csv(summary_tbl, file.path(out_dir, "summary.csv"), row.names = FALSE)

fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "f"))
md <- c(
  "# Pipeline summary", "",
  sprintf("Mode: `%s` | seed: %d | generated: %s", config$mode, config$seed, format(Sys.time(), "%Y-%m-%d %H:%M")), "",
  "## Cohorts", "",
  sprintf("- Discovery: n = %d (%s)", nrow(disc$x), paste(sprintf("%s = %d", levels(disc$y), table(disc$y)), collapse = ", ")),
  sprintf("- Validation: n = %d (%s)", nrow(vali$x), paste(sprintf("%s = %d", levels(vali$y), table(vali$y)), collapse = ", ")),
  sprintf("- CpGs shared between cohorts: %d", length(shared)), "",
  "## Stage 1: tri-filter (discovery cohort only)", "",
  sprintf("Parameters: variance keep fraction %s | significance %s, %s p < %s | |delta-beta| >= %s | top %d",
          tf$params$var_keep_frac, tf$params$sig_method,
          if (identical(tf$params$p_adjust, "none")) "unadjusted" else paste0(tf$params$p_adjust, "-adjusted"),
          tf$params$p_cutoff, tf$params$delta_beta_min, tf$params$top_k), "",
  "| Stage | CpGs remaining |", "|---|---|",
  sprintf("| Input | %d |", tf$counts[["input"]]),
  sprintf("| After variance filter | %d |", tf$counts[["after_variance"]]),
  sprintf("| After significance filter | %d |", tf$counts[["after_significance"]]),
  sprintf("| After effect-size filter | %d |", tf$counts[["after_effect"]]),
  sprintf("| Candidate pool (top k) | %d |", tf$counts[["final"]]),
  if (!is.null(truth)) c("", sprintf("Planted signal CpGs recovered in the candidate pool: %d of %d", sum(pool %in% truth), length(truth))) else NULL,
  "",
  "## Stage 2: permutation null (random panels from the candidate pool)", "",
  sprintf("%d random panels of %d CpGs, Random Forest with %dx%d repeated CV, scored on validation.",
          nrow(null), config$permutation$panel_size, config$permutation$cv_number, config$permutation$cv_repeats), "",
  "| Metric | Mean | SD | 95th percentile | Max |", "|---|---|---|---|---|",
  sprintf("| Validation AUC | %s | %s | %s | %s |", fmt(mean(null$valid_auc, na.rm = TRUE)), fmt(sd(null$valid_auc, na.rm = TRUE)),
          fmt(quantile(null$valid_auc, 0.95, na.rm = TRUE)), fmt(max(null$valid_auc, na.rm = TRUE))),
  sprintf("| Validation accuracy | %s | %s | %s | %s |", fmt(mean(null$valid_accuracy)), fmt(sd(null$valid_accuracy)),
          fmt(quantile(null$valid_accuracy, 0.95)), fmt(max(null$valid_accuracy))), "",
  "## Stages 3-4: RFE-selected panels on the validation cohort", "",
  "| Learner | CpGs | CV ROC | Val. AUC | Val. accuracy | Bal. accuracy | Sens. | Spec. | Empirical p (AUC vs null) |",
  "|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %d | %s | %s | %s | %s | %s | %s | %s |", summary_tbl$learner, summary_tbl$n_selected,
          fmt(summary_tbl$cv_roc), fmt(summary_tbl$valid_auc), fmt(summary_tbl$valid_accuracy),
          fmt(summary_tbl$valid_balanced_accuracy), fmt(summary_tbl$valid_sensitivity),
          fmt(summary_tbl$valid_specificity), fmt(summary_tbl$empirical_p_auc)), "",
  if (!is.null(truth)) sprintf("Planted signal CpGs in each selected panel: %s",
                               paste(sprintf("%s %d/%d", summary_tbl$learner, summary_tbl$planted_signal_recovered, summary_tbl$n_selected), collapse = "; ")) else NULL,
  "",
  sprintf("RFE profile (cross-validated ROC by panel size; chosen = smallest size within %s%% of the best):",
          config$rfe$size_tolerance %||% 3), "",
  vapply(rfe_res, function(r) {
    pr <- r$profile
    sprintf("- %s: %s (chosen %d)", r$learner,
            paste(sprintf("%d: %s", pr$Variables, fmt(pr$ROC)), collapse = " | "), length(r$selected))
  }, character(1)),
  "",
  "Read the columns with care. CV ROC is from the refit of an already-selected panel on the same",
  "discovery data, so it is optimistically biased; the validation columns are the honest estimate.",
  "The empirical p-value is the fraction of random same-pool panels whose validation AUC",
  "matched or exceeded the selected panel's (with +1 correction). A value near 0.5 means the",
  "selected panel is typical of random panels from the pool, i.e. RFE added nothing beyond",
  "the tri-filter; a small value means the selected panel beats chance draws from the same pool.",
  "Random panels have a fixed size, so a much larger selected panel has a size advantage."
)
writeLines(md[!is.na(md)], file.path(out_dir, "summary.md"))
saveRDS(list(config = config, trifilter = tf, permutation_null = null, rfe = rfe_res, summary = summary_tbl),
        file.path(out_dir, "pipeline_objects.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\n", paste(md, collapse = "\n"), "\n", sep = "")
log_msg("done: results in %s", out_dir)
