# 03_rfe_benchmark.R -- recursive feature elimination with three learners
#
# For each learner (Random Forest, linear SVM, radial SVM), caret::rfe evaluates
# nested panels of decreasing size inside repeated k-fold CV on the discovery
# cohort and picks the size that maximises cross-validated ROC. The selected
# panel is then refit with repeated CV and scored once on the validation cohort.
#
# Why this shape: RFE's resampling is *external* to the feature ranking, so the
# ranking is re-derived inside every fold. Ranking once on the full discovery
# set and then cross-validating the classifier would leak the selection into
# the CV estimate (Ambroise & McLachlan, 2002).

rfe_benchmark <- function(x_train, y_train, x_valid, y_valid,
                          learners = c("rf", "svmLinear", "svmRadial"),
                          sizes = 5:30, cv_number = 5, cv_repeats = 10,
                          ntree = 500, seed = 123, inner_cv_number = 3,
                          size_tolerance = 3) {
  sizes <- sort(unique(sizes[sizes < ncol(x_train)]))
  if (!length(sizes)) stop("No RFE subset sizes are smaller than the candidate pool (", ncol(x_train), " CpGs)")

  # Parsimony rule: take the SMALLEST panel whose CV ROC is within `size_tolerance`
  # percent of the best, rather than the single best size. rfe() always evaluates
  # the full pool too, and when the pool holds many weak true effects the curve
  # keeps rising, so "best" is simply "everything". A biomarker panel has to be
  # assayable, and the tolerance makes that trade-off explicit (0 = best size).
  pick_size <- function(x, metric, maximize) {
    caret::pickSizeTolerance(x, metric, tol = size_tolerance, maximize = maximize)
  }

  lapply(stats::setNames(learners, learners), function(learner) {
    log_msg("  RFE [%s]: subset sizes %s, %dx%d repeated CV, smallest size within %s%% of best ROC",
            learner, paste(range(sizes), collapse = "-"), cv_number, cv_repeats, size_tolerance)
    t0 <- Sys.time()
    set.seed(seed)

    if (learner == "rf") {
      funcs <- caret::rfFuncs
      funcs$summary    <- caret::twoClassSummary
      funcs$selectSize <- pick_size
      rctrl <- caret::rfeControl(functions = funcs, method = "repeatedcv",
                                 number = cv_number, repeats = cv_repeats, verbose = FALSE)
      rfe_fit <- caret::rfe(x = x_train, y = y_train, sizes = sizes,
                            rfeControl = rctrl, metric = "ROC", ntree = ntree)
    } else {
      funcs <- caret::caretFuncs
      funcs$summary    <- caret::twoClassSummary
      funcs$selectSize <- pick_size
      # rfe() keeps `metric` for itself, so the inner train() would otherwise fall
      # back to Accuracy (and warn on every fit); pin ROC explicitly.
      funcs$fit <- function(x, y, first, last, ...) caret::train(x, y, metric = "ROC", ...)
      rctrl <- caret::rfeControl(functions = funcs, method = "repeatedcv",
                                 number = cv_number, repeats = cv_repeats, verbose = FALSE)
      inner <- caret::trainControl(method = "cv", number = inner_cv_number,
                                   classProbs = TRUE, summaryFunction = caret::twoClassSummary)
      rfe_fit <- caret::rfe(x = x_train, y = y_train, sizes = sizes,
                            rfeControl = rctrl, metric = "ROC",
                            method = learner, trControl = inner, tuneLength = 3)
    }
    selected <- caret::predictors(rfe_fit)

    # Refit the chosen panel with repeated CV, then score it ONCE on validation.
    ctrl <- caret::trainControl(method = "repeatedcv", number = cv_number, repeats = cv_repeats,
                                classProbs = TRUE, summaryFunction = caret::twoClassSummary)
    set.seed(seed)
    final <- if (learner == "rf") {
      caret::train(x = x_train[, selected, drop = FALSE], y = y_train, method = "rf",
                   metric = "ROC", trControl = ctrl, ntree = ntree)
    } else {
      caret::train(x = x_train[, selected, drop = FALSE], y = y_train, method = learner,
                   metric = "ROC", trControl = ctrl, tuneLength = 5)
    }
    ev <- evaluate_model(final, x_valid[, selected, drop = FALSE], y_valid)
    log_msg("  RFE [%s]: %d CpGs selected | CV ROC %.3f | validation AUC %.3f, accuracy %.3f (%.1f min)",
            learner, length(selected), max(final$results$ROC), ev$auc, ev$accuracy,
            as.numeric(difftime(Sys.time(), t0, units = "mins")))

    list(learner = learner, rfe = rfe_fit, profile = rfe_fit$results,
         selected = selected, model = final,
         cv_roc = max(final$results$ROC), validation = ev)
  })
}
