# Unit tests for the pure functions. The end-to-end behaviour is covered by
# `Rscript scripts/run_pipeline.R --smoke` (see .github/workflows/smoke-test.yml).

test_that("make_outcome puts the positive class first and rejects stray labels", {
  y <- make_outcome(c("N", "P", "P", "N"), positive = "P", negative = "N")
  expect_identical(levels(y), c("P", "N"))
  expect_error(make_outcome(c("P", "N", "?"), "P", "N"), "unexpected labels")
})

test_that("beta_to_m is the log2 logit and never returns infinities", {
  expect_equal(beta_to_m(0.5), 0)
  expect_equal(beta_to_m(0.8), log2(4))
  expect_true(all(is.finite(beta_to_m(c(0, 1, 1e-9)))))
})

test_that("row_variance matches apply(var)", {
  set.seed(1)
  m <- matrix(rnorm(200), 20, 10)
  expect_equal(unname(row_variance(m)), apply(m, 1, var))
})

test_that("vectorised Welch test matches t.test row by row", {
  set.seed(2)
  m <- matrix(rnorm(30 * 12), 30, 12)
  is_pos <- rep(c(TRUE, FALSE), each = 6)
  p_vec <- welch_t_rows(m, is_pos)
  p_ref <- apply(m, 1, function(r) t.test(r[is_pos], r[!is_pos])$p.value)
  expect_equal(unname(p_vec), unname(p_ref), tolerance = 1e-10)
})

test_that("empirical_p applies the +1 correction and respects the alternative", {
  null <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  expect_equal(empirical_p(null, 0.85), (1 + 1) / (5 + 1))            # one null value >= 0.85
  expect_equal(empirical_p(null, 0.95), 1 / 6)                        # none, floor at 1/(n+1)
  expect_equal(empirical_p(null, 0.55, "less"), 2 / 6)                # one null value <= 0.55
  expect_equal(empirical_p(c(null, NA, Inf), 0.95), 1 / 6)            # non-finite values are dropped
})

test_that("tri-filter uses discovery data only and respects every threshold", {
  set.seed(3)
  n <- 60; p <- 300
  x <- matrix(runif(n * p, 0.2, 0.7), n, p, dimnames = list(sprintf("s%02d", 1:n), sprintf("cg%03d", 1:p)))
  y <- make_outcome(rep(c("P", "N"), each = n / 2), "P", "N")
  x[y == "P", 1:10] <- x[y == "P", 1:10] + 0.25                      # plant 10 strong CpGs
  x <- pmin(x, 0.999)
  res <- trifilter_select(x, y, var_keep_frac = 0.5, sig_method = "welch", p_cutoff = 0.05,
                          p_adjust = "BH", delta_beta_min = 0.05, top_k = 8)
  expect_equal(unname(res$counts[["after_variance"]]), 150)
  expect_length(res$features, 8)
  expect_true(all(res$features %in% sprintf("cg%03d", 1:10)))
  tab <- res$table
  expect_true(all(tab$p_adj[tab$selected] < 0.05))
  expect_true(all(abs(tab$delta_beta[tab$selected]) >= 0.05))
  expect_true(all(tab$cpg[tab$selected] %in% head(tab$cpg[order(tab$p_value)], 8)))
  expect_error(trifilter_select(x[, 1:5], y, sig_method = "welch", top_k = 80), NA)  # pool smaller than top_k is fine
})

test_that("tri-filter treats NA statistics as failing the filter", {
  set.seed(4)
  n <- 30
  x <- matrix(runif(n * 20, 0.3, 0.7), n, 20, dimnames = list(sprintf("s%02d", 1:n), sprintf("cg%03d", 1:20)))
  x[, 1] <- 0.5                                                       # zero-variance probe -> NA t statistic
  y <- make_outcome(rep(c("P", "N"), each = n / 2), "P", "N")
  res <- trifilter_select(x, y, var_keep_frac = 1, sig_method = "welch", p_adjust = "none",
                          p_cutoff = 0.5, delta_beta_min = 0, top_k = 20)
  expect_false("cg001" %in% res$features)
  expect_false(any(is.na(res$table$pass_significance)))
  expect_false(any(is.na(res$features)))
})

test_that("synthetic generator writes aligned cohorts with the planted signal", {
  out <- file.path(tempdir(), "synth-test")
  truth <- make_synthetic_methylation(out, n_discovery = 30, n_validation = 20, n_cpg = 500,
                                      n_signal = 10, delta_range = c(0.1, 0.2), seed = 7)
  disc <- load_cohort(file.path(out, "discovery_beta.csv"), file.path(out, "discovery_labels.csv"),
                      "sample_id", "status", "P", "N")
  expect_equal(dim(disc$x), c(30, 500))
  expect_identical(levels(disc$y), c("P", "N"))
  expect_setequal(readLines(file.path(out, "truth_signal_cpgs.txt")), truth$signal)
  delta <- colMeans(disc$x[disc$y == "P", truth$signal]) - colMeans(disc$x[disc$y == "N", truth$signal])
  expect_equal(sign(delta), sign(truth$delta[truth$signal]))
  expect_true(all(disc$x > 0 & disc$x < 1))
})
