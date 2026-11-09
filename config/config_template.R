# config_template.R -- copy to config/my_study.R, edit, then run:
#   Rscript scripts/run_pipeline.R --config config/my_study.R
#
# Input format
#   *_beta.csv    CpG x sample matrix; first column = CpG ID, header = sample IDs
#   *_labels.csv  one row per sample with at least `sample_col` and `outcome_col`
#
# Feature selection runs on the discovery cohort only. The validation cohort
# is used exactly once per model, for scoring.

config <- list(
  out_dir = "results/my_study",

  data = list(
    synthetic         = FALSE,
    discovery_beta    = "data/discovery_beta.csv",
    discovery_labels  = "data/discovery_labels.csv",
    validation_beta   = "data/validation_beta.csv",
    validation_labels = "data/validation_labels.csv",
    truth_signal      = NULL,          # optional: file of known signal CpGs (simulations only)
    sample_col  = "sample_id",
    outcome_col = "status",
    positive    = "P",                 # the event class; ROC/Sens/Spec refer to this class
    negative    = "N"
  ),

  # Stage 1 -- tri-filter. These defaults are placeholders: set them to the
  # thresholds used in your study and report them in the methods.
  #
  # An unadjusted p-value is a deliberately liberal pre-filter. It is defensible
  # here only because selection never touches the validation cohort and the
  # permutation null measures what the resulting pool gives by chance. Set
  # p_adjust = "BH" (and a cutoff such as 0.05) if the pool itself must be an
  # FDR-controlled list; at typical EWAS sample sizes that pool may be empty.
  trifilter = list(
    var_keep_frac  = 0.25,             # keep the top 25% most variable CpGs
    sig_method     = "limma",          # "limma" (moderated t on M-values) or "welch" (t-test on betas)
    p_cutoff       = 0.01,
    p_adjust       = "none",           # any p.adjust() method, or "none"
    delta_beta_min = 0.05,             # minimum |mean beta difference| between classes
    top_k          = 80                # size of the candidate pool handed to RFE
  ),

  # Stage 2 -- permutation null. 1000 x (10x10 CV) is the full design; it is
  # CPU-heavy, so raise n_cores or lower cv_repeats for a first pass.
  permutation = list(
    panel_size = 10,
    n_perm     = 1000,
    cv_number  = 10,
    cv_repeats = 10,
    ntree      = 500
  ),

  # Stage 3 -- recursive feature elimination. rfe() always evaluates the full
  # pool as well as `sizes`. With many weak true effects, CV ROC keeps rising
  # with panel size, so the chosen size is an assay-cost decision: the smallest
  # panel within `size_tolerance` percent of the best CV ROC (0 = best size).
  rfe = list(
    learners       = c("rf", "svmLinear", "svmRadial"),
    sizes          = c(5, 10, 15, 20, 25, 30, 40, 50),
    cv_number      = 5,
    cv_repeats     = 10,
    ntree          = 500,
    size_tolerance = 3
  ),

  seed    = 123,
  n_cores = 1
)
