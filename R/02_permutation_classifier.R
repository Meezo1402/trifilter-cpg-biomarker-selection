# 02_permutation_classifier.R -- empirical null distribution from random CpG panels
#
# Draws n_perm random panels of panel_size CpGs from the post-filter candidate
# pool, trains a Random Forest with repeated k-fold CV on the discovery cohort,
# and scores every panel on the held-out validation cohort.
#
# The resulting distribution is the null against which any *selected* panel is
# judged. A panel chosen by recursive feature elimination has to outperform
# random panels drawn from the same pool; if it does not, the selection step
# added nothing beyond what the tri-filter already achieved. This separates
# "the filters found signal" from "RFE found a better panel than chance".

permutation_null <- function(x_train, y_train, x_valid, y_valid, feature_pool,
                             panel_size = 10, n_perm = 1000,
                             cv_number = 10, cv_repeats = 10,
                             ntree = 500, seed = 123, n_cores = 1) {
  if (length(feature_pool) < panel_size) {
    stop("Candidate pool (", length(feature_pool), ") is smaller than panel_size (", panel_size, ")")
  }

  one_perm <- function(i) {
    set.seed(seed + i)   # per-permutation seed: reproducible regardless of core scheduling
    panel <- sample(feature_pool, panel_size)
    ctrl  <- caret::trainControl(method = "repeatedcv", number = cv_number, repeats = cv_repeats,
                                 classProbs = TRUE, summaryFunction = caret::twoClassSummary)
    fit <- caret::train(x = x_train[, panel, drop = FALSE], y = y_train, method = "rf",
                        metric = "ROC", trControl = ctrl, ntree = ntree,
                        tuneGrid = data.frame(mtry = max(1L, floor(sqrt(panel_size)))))
    ev <- evaluate_model(fit, x_valid[, panel, drop = FALSE], y_valid)
    data.frame(perm_id = i, panel = collapse_panel(panel),
               cv_roc = max(fit$results$ROC),
               valid_auc = ev$auc,
               valid_accuracy = ev$accuracy,
               valid_balanced_accuracy = ev$balanced_accuracy,
               valid_sensitivity = ev$sensitivity,
               valid_specificity = ev$specificity,
               stringsAsFactors = FALSE)
  }

  log_msg("  permutation null: %d random panels of %d CpGs, RF with %dx%d repeated CV, %d core(s)",
          n_perm, panel_size, cv_number, cv_repeats, n_cores)
  t0 <- Sys.time()
  res <- if (n_cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(n_perm), one_perm, mc.cores = n_cores)
  } else {
    lapply(seq_len(n_perm), one_perm)
  }
  failed <- vapply(res, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(sum(failed), " permutation(s) failed; first error: ", as.character(res[[which(failed)[1]]]))
  res <- do.call(rbind, res)
  log_msg("  permutation null finished in %.1f min; validation AUC mean %.3f, 95th percentile %.3f",
          as.numeric(difftime(Sys.time(), t0, units = "mins")),
          mean(res$valid_auc, na.rm = TRUE), stats::quantile(res$valid_auc, 0.95, na.rm = TRUE))
  res
}
