# 01_trifilter.R -- three-stage CpG pre-selection on the DISCOVERY cohort only
#
#   Stage 1  Variance filter      keep the most variable CpGs; near-constant probes
#                                 carry no discriminative information and inflate
#                                 the multiple-testing burden
#   Stage 2  Significance filter  differential methylation between classes
#                                 (limma moderated t on M-values, or vectorised
#                                 Welch t on beta values), multiplicity-adjusted
#   Stage 3  Effect-size filter   minimum absolute mean beta difference (delta-beta),
#                                 so that statistically significant but biologically
#                                 trivial shifts are discarded
#   Final    Ranking              survivors ranked by p-value (ties by |delta-beta|);
#                                 the top_k CpGs form the candidate pool for RFE
#
# The validation cohort is never seen by any stage. Feature selection performed
# on pooled data is the single most common source of optimistic bias in
# biomarker studies; keeping it discovery-only is deliberate.

trifilter_select <- function(x, y,
                             var_keep_frac  = 0.25,
                             sig_method     = c("limma", "welch"),
                             p_cutoff       = 0.05,
                             p_adjust       = "BH",
                             delta_beta_min = 0.05,
                             top_k          = 80) {
  sig_method <- match.arg(sig_method)
  stopifnot(is.matrix(x), is.factor(y), nlevels(y) == 2, nrow(x) == length(y))
  positive <- levels(y)[1]
  negative <- levels(y)[2]
  beta <- t(x)                       # CpG x sample for row-wise statistics
  n_input <- nrow(beta)

  # ---- Stage 1: variance ---------------------------------------------------
  v <- row_variance(beta)
  keep_n <- max(1L, floor(var_keep_frac * n_input))
  stage1 <- names(sort(v, decreasing = TRUE))[seq_len(keep_n)]
  beta1 <- beta[stage1, , drop = FALSE]

  # ---- Stage 2: significance ----------------------------------------------
  is_pos   <- y == positive
  mean_pos <- rowMeans(beta1[, is_pos, drop = FALSE])
  mean_neg <- rowMeans(beta1[, !is_pos, drop = FALSE])
  delta_beta <- mean_pos - mean_neg

  if (sig_method == "limma" && !requireNamespace("limma", quietly = TRUE)) {
    log_msg("  limma not installed; falling back to the Welch t-test")
    sig_method <- "welch"
  }
  if (sig_method == "limma") {
    m      <- beta_to_m(beta1)
    group  <- factor(ifelse(is_pos, positive, negative), levels = c(negative, positive))
    design <- model.matrix(~ group)                      # coefficient 2 = positive vs negative
    fit    <- limma::eBayes(limma::lmFit(m, design))
    tt     <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
    idx    <- match(rownames(beta1), rownames(tt))
    p_value  <- tt$P.Value[idx]
    log2fc_m <- tt$logFC[idx]
  } else {
    p_value  <- welch_t_rows(beta1, is_pos)
    log2fc_m <- NA_real_
  }
  p_adj <- if (identical(p_adjust, "none")) p_value else p.adjust(p_value, method = p_adjust)

  tab <- data.frame(cpg = rownames(beta1), variance = unname(v[stage1]),
                    mean_pos = unname(mean_pos), mean_neg = unname(mean_neg),
                    delta_beta = unname(delta_beta), log2fc_m = log2fc_m,
                    p_value = p_value, p_adj = p_adj,
                    stringsAsFactors = FALSE)
  tab$pass_significance <- !is.na(tab$p_adj) & tab$p_adj < p_cutoff   # NA p (e.g. zero-variance probe) never passes

  # ---- Stage 3: effect size -----------------------------------------------
  tab$pass_effect <- !is.na(tab$delta_beta) & abs(tab$delta_beta) >= delta_beta_min
  survivors <- tab[tab$pass_significance & tab$pass_effect, ]

  # ---- Ranking --------------------------------------------------------------
  survivors <- survivors[order(survivors$p_value, -abs(survivors$delta_beta)), ]
  selected  <- head(survivors$cpg, top_k)
  tab$selected <- tab$cpg %in% selected

  counts <- c(input = n_input, after_variance = length(stage1),
              after_significance = sum(tab$pass_significance),
              after_effect = nrow(survivors), final = length(selected))
  log_msg("  tri-filter: %d -> variance %d -> significance %d -> effect size %d -> top %d",
          counts[["input"]], counts[["after_variance"]], counts[["after_significance"]],
          counts[["after_effect"]], counts[["final"]])

  list(features = selected, table = tab, counts = counts,
       params = list(var_keep_frac = var_keep_frac, sig_method = sig_method,
                     p_cutoff = p_cutoff, p_adjust = p_adjust,
                     delta_beta_min = delta_beta_min, top_k = top_k))
}

#' Vectorised Welch two-sample t-test across all rows at once.
welch_t_rows <- function(m, is_pos) {
  a <- m[, is_pos, drop = FALSE]
  b <- m[, !is_pos, drop = FALSE]
  na <- ncol(a); nb <- ncol(b)
  ma <- rowMeans(a); mb <- rowMeans(b)
  va <- row_variance(a); vb <- row_variance(b)
  se2 <- va / na + vb / nb
  t_stat <- (ma - mb) / sqrt(se2)
  df <- se2^2 / ((va / na)^2 / (na - 1) + (vb / nb)^2 / (nb - 1))
  2 * stats::pt(-abs(t_stat), df)
}
