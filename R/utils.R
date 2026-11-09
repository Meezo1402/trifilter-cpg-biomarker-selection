# utils.R -- shared helpers for the tri-filter biomarker pipeline
# Author: Mazen Istanbouli

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
}

#' Stop early with an actionable message if required packages are missing.
check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         "\nInstall with: install.packages(c(",
         paste(sprintf('"%s"', missing), collapse = ", "), "))",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Row-wise variance, vectorised (no matrixStats dependency).
row_variance <- function(m) {
  rowSums((m - rowMeans(m))^2) / (ncol(m) - 1)
}

#' Beta values -> M-values, clipped to avoid infinite logits at 0 and 1.
beta_to_m <- function(beta, eps = 1e-3) {
  beta <- pmin(pmax(beta, eps), 1 - eps)
  log2(beta / (1 - beta))
}

#' Build the binary outcome factor with the POSITIVE class as the first level.
#' caret::twoClassSummary treats the first factor level as the event of
#' interest, so this ordering is what makes ROC/Sens/Spec refer to the
#' positive class everywhere downstream.
make_outcome <- function(y, positive, negative) {
  y <- as.character(y)
  bad <- setdiff(unique(y), c(positive, negative))
  if (length(bad)) {
    stop("Outcome contains unexpected labels: ", paste(bad, collapse = ", "),
         " (expected '", positive, "' or '", negative, "')")
  }
  factor(y, levels = c(positive, negative))
}

#' Load one cohort: a CpG x sample beta matrix (CSV, first column = CpG ID)
#' and a label table. Samples are intersected and aligned; the matrix is
#' returned transposed (samples x CpGs) for modelling.
load_cohort <- function(beta_path, labels_path, sample_col, outcome_col,
                        positive, negative) {
  beta <- as.matrix(read.csv(beta_path, row.names = 1, check.names = FALSE))
  labs <- read.csv(labels_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c(sample_col, outcome_col) %in% colnames(labs))) {
    stop("Label file must contain the columns '", sample_col, "' and '", outcome_col, "'")
  }
  common <- intersect(colnames(beta), labs[[sample_col]])
  if (!length(common)) stop("No sample IDs shared between beta matrix columns and the label file")
  dropped <- length(union(colnames(beta), labs[[sample_col]])) - length(common)
  if (dropped > 0) log_msg("  note: %d sample(s) dropped (present in only one of beta/labels)", dropped)
  beta <- beta[, common, drop = FALSE]
  labs <- labs[match(common, labs[[sample_col]]), , drop = FALSE]
  list(x = t(beta),
       y = make_outcome(labs[[outcome_col]], positive, negative),
       sample_ids = common)
}

#' Collapse a feature panel to one string for CSV output.
collapse_panel <- function(v) paste(v, collapse = ";")

#' Empirical p-value with the +1 correction (Phipson & Smyth, 2010):
#' the probability that a random panel performs at least as well as the
#' observed panel.
empirical_p <- function(null_values, observed, alternative = c("greater", "less")) {
  alternative <- match.arg(alternative)
  null_values <- null_values[is.finite(null_values)]
  hits <- if (alternative == "greater") sum(null_values >= observed) else sum(null_values <= observed)
  (hits + 1) / (length(null_values) + 1)
}
