# 04_evaluate.R -- score a fitted caret model on the held-out validation cohort
#
# Reports threshold-free (AUC) and threshold-dependent (accuracy, balanced
# accuracy, sensitivity, specificity, kappa) metrics. Balanced accuracy is
# reported alongside raw accuracy because class imbalance in validation cohorts
# can make raw accuracy look good for a classifier that predicts one class.
# If the validation cohort contains a single class, AUC is NA rather than a
# misleading number.

evaluate_model <- function(model, x_valid, y_valid) {
  positive <- levels(y_valid)[1]
  negative <- levels(y_valid)[2]
  pred_class <- predict(model, newdata = x_valid)
  pred_prob  <- predict(model, newdata = x_valid, type = "prob")[[positive]]
  cm <- caret::confusionMatrix(pred_class, y_valid, positive = positive)

  auc <- if (nlevels(droplevels(y_valid)) == 2) {
    as.numeric(pROC::auc(pROC::roc(response = y_valid, predictor = pred_prob,
                                   levels = c(negative, positive), direction = "<",
                                   quiet = TRUE)))
  } else {
    NA_real_
  }

  list(confusion = cm,
       auc = auc,
       accuracy = unname(cm$overall["Accuracy"]),
       balanced_accuracy = unname(cm$byClass["Balanced Accuracy"]),
       sensitivity = unname(cm$byClass["Sensitivity"]),
       specificity = unname(cm$byClass["Specificity"]),
       kappa = unname(cm$overall["Kappa"]),
       predictions = data.frame(sample = rownames(x_valid), observed = y_valid,
                                predicted = pred_class, prob_positive = pred_prob,
                                stringsAsFactors = FALSE))
}
