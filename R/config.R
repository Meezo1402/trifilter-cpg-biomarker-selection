# config.R -- built-in configurations for the synthetic demo and the CI smoke test
#
# Your own study is configured by copying config/config_template.R and passing
# it with --config; these two presets exist so the pipeline can be run with
# zero setup.

default_config <- function(root, mode = c("demo", "smoke")) {
  mode <- match.arg(mode)
  data_dir <- file.path(root, "data", "synthetic", mode)
  n_cores  <- max(1L, min(4L, parallel::detectCores() - 1L))

  base <- list(
    mode    = mode,
    out_dir = file.path(root, "results", mode),
    data = list(
      synthetic         = TRUE,                     # generated on first run if absent
      discovery_beta    = file.path(data_dir, "discovery_beta.csv"),
      discovery_labels  = file.path(data_dir, "discovery_labels.csv"),
      validation_beta   = file.path(data_dir, "validation_beta.csv"),
      validation_labels = file.path(data_dir, "validation_labels.csv"),
      truth_signal      = file.path(data_dir, "truth_signal_cpgs.txt"),
      sample_col = "sample_id", outcome_col = "status",
      positive = "P", negative = "N"
    ),
    seed = 123, n_cores = n_cores
  )

  if (mode == "demo") {
    # Modest signal on purpose (see synthetic_data.R): 100 planted CpGs with
    # |delta-beta| 0.02-0.08, of which only the stronger third or so reach
    # unadjusted p < 0.01 at n = 186. The candidate pool therefore mixes true and
    # false positives, which is the regime a liberal pre-filter + independent
    # validation design is meant for. Larger shifts make every panel score
    # AUC ~1 and the benchmark stops saying anything.
    base$synthetic   <- list(n_discovery = 186, n_validation = 90, n_cpg = 20000, n_signal = 100,
                             delta_range = c(0.02, 0.08), noise_sd = 0.5)
    base$trifilter   <- list(var_keep_frac = 0.25, sig_method = "limma", p_cutoff = 0.01,
                             p_adjust = "none", delta_beta_min = 0.05, top_k = 80)
    base$permutation <- list(panel_size = 10, n_perm = 200, cv_number = 5, cv_repeats = 3, ntree = 500)
    base$rfe         <- list(learners = c("rf", "svmLinear", "svmRadial"),
                             sizes = c(5, 10, 15, 20, 25, 30, 40, 50), cv_number = 5, cv_repeats = 3,
                             ntree = 500, size_tolerance = 3)
  } else {
    # Strong signal so the smoke test never fails for lack of survivors.
    base$synthetic   <- list(n_discovery = 60, n_validation = 30, n_cpg = 2000, n_signal = 20,
                             delta_range = c(0.08, 0.15), noise_sd = 0.45)
    base$trifilter   <- list(var_keep_frac = 0.5, sig_method = "welch", p_cutoff = 0.05,
                             p_adjust = "none", delta_beta_min = 0.03, top_k = 30)
    base$permutation <- list(panel_size = 5, n_perm = 10, cv_number = 3, cv_repeats = 1, ntree = 100)
    base$rfe         <- list(learners = c("rf", "svmRadial"),
                             sizes = c(5, 10, 15), cv_number = 3, cv_repeats = 1, ntree = 100,
                             size_tolerance = 3)
    base$n_cores <- 1L
  }
  base
}
