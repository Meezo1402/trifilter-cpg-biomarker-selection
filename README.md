# trifilter_cpg_biomarker_selection
R pipeline for compact CpG biomarker panels from Illumina 450K/EPIC methylation data: a variance, significance and effect-size tri-filter on a discovery cohort, RF/SVM recursive feature elimination with nested CV, a permutation null of random same-pool panels, and one-shot scoring on an independent validation cohort. Synthetic demo, unit tests, CI.
An R pipeline for selecting a compact CpG biomarker panel from Illumina
450K/EPIC DNA methylation data and validating it on an independent cohort,
with an empirical permutation null that tests whether the selection step
beats chance draws from the same candidate pool.

Developed for a discovery cohort of 186 veterans with an independent
validation sample of 90 (Illumina 450K). The study data are under consortium
agreements and are not included; the repository runs end to end on synthetic
cohorts with planted signal, so every stage can be inspected and checked for
signal recovery.

## Why this design

Biomarker-selection studies fail in predictable ways: features are chosen on
pooled data and then "validated" on part of it, classifier accuracy is
reported without asking what a random panel would score, and a single
cross-validated number stands in for generalisation. Each stage below exists
to close one of those gaps.

```
discovery cohort                              validation cohort
      |                                              |
      v                                              |
 [1] TRI-FILTER   variance -> significance -> effect size -> top-k pool
      |                                              |
      +---> [2] PERMUTATION NULL                     |
      |         n random panels from the pool,       |
      |         RF + repeated CV ----- score -------->|   distribution of "chance" performance
      |                                              |
      +---> [3] RFE  (RF, linear SVM, radial SVM)    |
                nested resampling; size = smallest    |
                within tolerance of the best CV ROC   |
                refit selected panel ----- score ---->|   AUC, accuracy, sens/spec
                                                     |
            [4] EMPIRICAL P  =  P(random panel >= selected panel)
```

| Stage | What it does | Why |
|---|---|---|
| 1. Tri-filter | Variance filter, then differential methylation (limma moderated *t* on M-values, or Welch *t*; adjusted or unadjusted *p*), then a minimum delta-beta; survivors ranked by *p* and the top *k* kept | Reduces 450K probes to a tractable candidate pool using **discovery data only**; the effect-size gate drops significant-but-trivial shifts. A liberal *p* cut-off is acceptable here precisely because the pool is never judged on the data that produced it |
| 2. Permutation null | Random panels of *m* CpGs drawn from the pool, each trained with repeated CV and scored on validation | Establishes what any panel from this pool achieves by chance, so a selected panel is judged against the pool it came from, not against 50% |
| 3. RFE benchmark | `caret::rfe` with RF, linear SVM, and radial SVM; ranking re-derived inside every resample; panel size = the smallest within 3% of the best cross-validated ROC | Avoids the selection-bias leak of ranking once then cross-validating (Ambroise & McLachlan, 2002); the parsimony rule makes the panel-size / assay-cost trade-off explicit instead of defaulting to the whole pool |
| 4. Validation | Each selected panel refit with repeated CV and scored **once** on the held-out cohort; empirical *p* against the null | Threshold-free (AUC) and threshold-dependent metrics, with balanced accuracy for imbalanced cohorts |

## Quick start

```bash
git clone https://github.com/meezo1402/trifilter-cpg-biomarker-selection.git
cd trifilter-cpg-biomarker-selection
Rscript -e 'install.packages(c("caret", "randomForest", "kernlab", "pROC", "e1071"))'
# optional, for the limma significance filter (falls back to Welch t if absent):
Rscript -e 'BiocManager::install("limma")'

Rscript scripts/run_pipeline.R --demo     # 20,000 CpGs, n = 186 / 90, ~5 min on one core
Rscript scripts/run_pipeline.R --smoke    # tiny cohort, ~10 s; what CI runs
Rscript tests/testthat.R                  # unit tests (needs testthat)
```

Results land in `results/<mode>/`: a per-CpG table with every filter decision,
the permutation null, each learner's RFE profile and selected panel,
per-sample validation predictions, `summary.md`, and `sessionInfo.txt`.

## Running on your own cohorts

```bash
cp config/config_template.R config/my_study.R   # edit paths and thresholds
Rscript scripts/run_pipeline.R --config config/my_study.R
```

Inputs are a CpG x sample beta matrix (CSV) and a label table per cohort; see
`data/README.md`. The filter thresholds in the template are placeholders and
should be set to, and reported as, the values used in your study. The full
permutation design (1000 panels x 10x10 repeated CV) is CPU-heavy; set
`n_cores` in the config to parallelise it.

## Example output (synthetic demo)

Numbers below are from `--demo` on synthetic cohorts (seed 123) and say
nothing about the study data. The generator plants 100 CpGs with
|delta-beta| 0.02-0.08 among 20,000; at n = 186 only the stronger third reach
unadjusted *p* < 0.01, so the pool is a realistic mix of true and false
positives and the downstream stages have something to do.

**Stage 1, tri-filter** (variance top 25% | limma, unadjusted *p* < 0.01 | |delta-beta| >= 0.05 | top 80)

| Stage | CpGs remaining |
|---|---|
| Input | 20,000 |
| After variance filter | 5,000 |
| After significance filter | 95 |
| After effect-size filter | 59 |
| Candidate pool | 59 (36 planted, 23 false positives) |

**Stage 2, permutation null** (200 random 10-CpG panels, RF, 5x3 CV, scored on validation)

| Metric | Mean | SD | 95th pct | Max |
|---|---|---|---|---|
| Validation AUC | 0.804 | 0.056 | 0.897 | 0.924 |
| Validation accuracy | 0.718 | 0.054 | 0.812 | 0.856 |

**Stages 3-4, selected panels on the held-out cohort**

| Learner | CpGs | planted | CV ROC (optimistic) | Val. AUC | Bal. accuracy | Sens. | Spec. | Empirical *p* |
|---|---|---|---|---|---|---|---|---|
| RF | 30 | 25 | 0.983 | 0.966 | 0.879 | 0.851 | 0.907 | 0.005 |
| linear SVM | 40 | 31 | 0.988 | 0.969 | 0.911 | 0.915 | 0.907 | 0.005 |
| radial SVM | 30 | 25 | 0.994 | 0.951 | 0.910 | 0.936 | 0.884 | 0.005 |

None of the 200 random same-pool panels reached the selected panels' validation
AUC (*p* = 1/201). The gap between CV ROC and validation AUC is the selection
optimism the design notes warn about. RFE's size-vs-ROC profile is in
`summary.md`; with 36 weak-to-moderate true CpGs in the pool it rises all the
way to the full pool (RF 0.985 at 59 CpGs), which is why a parsimony
tolerance rather than "best size" decides the panel. Whole run: about 5 min on
one core.

## Repository layout

```
R/
  01_trifilter.R               three-stage CpG pre-selection (discovery only)
  02_permutation_classifier.R  random-panel null distribution
  03_rfe_benchmark.R           RFE with RF / linear SVM / radial SVM
  04_evaluate.R                validation metrics (AUC, accuracy, balanced accuracy, sens/spec)
  synthetic_data.R             array-like synthetic cohorts with planted signal
  config.R                     built-in demo and smoke presets
  utils.R                      loading, outcome factor, empirical p, helpers
scripts/
  run_pipeline.R               end-to-end driver (--demo | --smoke | --config FILE)
  make_synthetic_data.R        standalone benchmark generator
tests/                         unit tests (Welch t vs t.test, empirical p, filter thresholds, generator)
config/config_template.R       template for your own study
data/README.md                 data availability and input format
.github/workflows/             unit tests + smoke run on every push
```

## Design notes

- **The outcome factor puts the positive class first.** `caret::twoClassSummary`
  treats the first level as the event, so ROC, sensitivity, and specificity
  refer to the positive class consistently across every stage.
- **Seeds are per permutation** (`seed + i`), so the null is reproducible
  whether it runs on one core or eight.
- **Validation is scored once per model.** Nothing from the validation cohort
  feeds back into filtering, ranking, or tuning.
- **The refit's CV ROC is optimistic and is labelled as such.** Once a panel
  has been selected on the discovery data, cross-validating the refit on that
  same data inherits the selection; the validation columns are the honest
  estimate, and the summary says so next to the table.
- **Panel size is a stated trade-off, not a hidden one.** With many weak true
  effects the CV curve keeps rising with size, so "best" would always be the
  whole pool. `size_tolerance` (default 3%) picks the smallest panel within
  that margin of the best, and the full size-vs-ROC profile is written out so
  a reader can disagree with the choice.
- **AUC is `NA`, not a number, if a validation cohort has one class.** Accuracy
  on a single-class cohort is reported but is not evidence of discrimination.
- **The null is fixed-size by design.** Random panels have `panel_size` CpGs
  (10 in the original study). A selected panel that is much larger has a size
  advantage in that comparison, so its empirical *p* is anti-conservative;
  compare like with like when the sizes differ a lot.
- **NA statistics fail the filter.** A zero-variance probe gives an NA *t*
  statistic; it is dropped rather than silently carried into the pool.
- **No data paths are hard-coded.** Cohorts, labels, thresholds, and compute
  are all config; the code never needs editing to run a new study.

## Dependencies

R >= 4.1; `caret`, `randomForest`, `kernlab`, `pROC`, `e1071`; optionally
`limma` (Bioconductor). Parallelism uses base `parallel` (fork-based; on
Windows the permutation stage runs serially).

## Provenance

This repository is a parameterised reimplementation of the analysis scripts
written for the original study (R, 2021). The study-specific data loaders,
cohort labels, and hard-coded model formulas were removed, the filter
thresholds were made configuration, and the leakage controls, seeding,
parsimony rule, and tests described above were added when the code was
prepared for release. The effect-size filter and the limma path assume
beta values; for other modalities use `sig_method = "welch"` on an
appropriately transformed matrix and set `delta_beta_min` on that scale.

## Citing

> Istanbouli, M. (2026). *trifilter-cpg-biomarker-selection: permutational
> tri-filter feature selection for methylation biomarkers.* GitHub repository,
> https://github.com/meezo1402/trifilter-cpg-biomarker-selection

## License

MIT. See `LICENSE`.
