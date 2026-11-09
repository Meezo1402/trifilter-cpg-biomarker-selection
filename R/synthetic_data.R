# synthetic_data.R -- simulated methylation cohorts for demonstration and testing
#
# The cohorts this pipeline was developed on are governed by consortium data
# agreements and cannot be redistributed. This generator produces array-like
# beta matrices with a known set of planted signal CpGs so that the full
# pipeline can be run, inspected, and checked for signal recovery by anyone.
#
# Design: baseline methylation is bimodal (most CpGs near 0 or 1, as on real
# arrays), sample-level noise is added on the logit scale, and n_signal CpGs
# receive a class-dependent beta shift in a random direction.
#
# The defaults are deliberately modest: 0.5% of CpGs carry a |delta-beta| of
# 0.02-0.08, and only the stronger ones are detectable at n = 186. Larger shifts
# make every panel score AUC ~1 and the benchmark stops saying anything.

make_synthetic_methylation <- function(out_dir,
                                       n_discovery = 186, n_validation = 90,
                                       n_cpg = 20000, n_signal = 100,
                                       delta_range = c(0.02, 0.08), noise_sd = 0.5,
                                       positive_frac = 0.5,
                                       positive = "P", negative = "N", seed = 42) {
  set.seed(seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cpg_ids <- sprintf("cg%08d", seq_len(n_cpg))

  base_beta  <- pmin(pmax(stats::rbeta(n_cpg, 0.5, 0.5), 0.02), 0.98)
  base_logit <- log(base_beta / (1 - base_beta))
  # Plant signal on intermediately methylated CpGs: a beta shift at a CpG sitting
  # near 0 or 1 is clipped away, and real differentially methylated positions
  # are overwhelmingly found in the intermediate range.
  eligible <- cpg_ids[base_beta > 0.2 & base_beta < 0.8]
  signal <- sample(eligible, n_signal)
  delta  <- stats::runif(n_signal, delta_range[1], delta_range[2]) *
            sample(c(-1, 1), n_signal, replace = TRUE)
  names(delta) <- signal

  simulate <- function(n, prefix) {
    y   <- sample(c(positive, negative), n, replace = TRUE, prob = c(positive_frac, 1 - positive_frac))
    ids <- sprintf("%s%03d", prefix, seq_len(n))
    logit <- matrix(base_logit, nrow = n_cpg, ncol = n) +
             matrix(stats::rnorm(n_cpg * n, 0, noise_sd), n_cpg, n)
    beta <- 1 / (1 + exp(-logit))
    dimnames(beta) <- list(cpg_ids, ids)
    is_pos <- y == positive
    beta[signal, is_pos] <- beta[signal, is_pos] + delta     # recycles delta down each column
    beta <- pmin(pmax(beta, 0.001), 0.999)
    dimnames(beta) <- list(cpg_ids, ids)                     # pmin/pmax strip attributes
    list(beta = beta, labels = data.frame(sample_id = ids, status = y, stringsAsFactors = FALSE))
  }

  disc <- simulate(n_discovery, "D")
  vali <- simulate(n_validation, "V")
  write.csv(disc$beta,   file.path(out_dir, "discovery_beta.csv"))
  write.csv(disc$labels, file.path(out_dir, "discovery_labels.csv"), row.names = FALSE)
  write.csv(vali$beta,   file.path(out_dir, "validation_beta.csv"))
  write.csv(vali$labels, file.path(out_dir, "validation_labels.csv"), row.names = FALSE)
  writeLines(signal, file.path(out_dir, "truth_signal_cpgs.txt"))
  log_msg("  synthetic data written to %s: %d CpGs, %d planted signal CpGs, discovery n=%d, validation n=%d",
          out_dir, n_cpg, n_signal, n_discovery, n_validation)
  invisible(list(signal = signal, delta = delta))
}
