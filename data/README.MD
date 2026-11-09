# Data

No study data are included in this repository.

The pipeline was developed on Illumina 450K DNA methylation data from a
discovery cohort of military veterans with an independent validation cohort.
Those data are governed by consortium data-use and publication agreements and
cannot be redistributed here.

Everything in this repository runs on **synthetic cohorts** that reproduce the
structure of array methylation data (bimodal baseline beta values, sample-level
noise, a planted set of differentially methylated CpGs). They are generated on
first run:

```
Rscript scripts/run_pipeline.R --demo          # writes data/synthetic/demo/
Rscript scripts/make_synthetic_data.R out_dir  # generate a custom benchmark
```

Because the planted signal is known, the demo reports how many true CpGs the
tri-filter recovers and how many end up in each selected panel. Real data have
no such ground truth; the permutation null is what stands in for it.

To run on your own cohorts, see `config/config_template.R`. Expected input:

| File | Format |
|---|---|
| `*_beta.csv` | CpG x sample matrix; first column CpG ID, header row = sample IDs |
| `*_labels.csv` | one row per sample; must contain the sample-ID and outcome columns named in the config |

`data/synthetic/` is git-ignored.
