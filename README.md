# CABNB

CABNB fits a compositionality-aware negative-binomial differential-abundance
model with empirical Bayes shrinkage, prevalence inference, and a group-blind
information filter.

## Installation

```r
devtools::install_github("noriakis/CABNB")
```

## Fit, inspect, and plot

```r
library(CABNB)

data(cabnb_example_counts)
data(cabnb_example_metadata)

fit <- cabnb_fit(
  counts = cabnb_example_counts,
  metadata = cabnb_example_metadata,
  formula = ~ group,
  coef = "grouptreatment"
)

fit
head(fit$results)
cabnb_discoveries(fit, type = "abundance")
cabnb_discoveries(fit, type = "information_filtered")
plot(fit)
```

The returned `fit$results` contains the abundance and prevalence estimates,
their p- and q-values, shrunken fold changes, and the information-filtered
test results. The main filtered decision columns are `information_eligible`,
`information_qvalue`, and `information_reject`.

## Run one MIDASim condition

The repository includes `scripts/midasim_one_condition.R`, which simulates one
two-group condition, fits CABNB, and calculates abundance, prevalence, and
combined-signal performance metrics. Run the following from the repository
root. The `MIDASim` and `PRROC` packages are required.

```r
library(CABNB)
source("scripts/midasim_one_condition.R")

result <- cabnb_run_midasim_one_condition(
  template = "throat",
  n_per_group = 40,
  n_taxa = 80,
  da_prop = 0.10,
  prevalence_prop = 0.05,
  effect_log_fc = log(2),
  prevalence_drop = 0.5,
  seed = 1,
  save_prefix = "results/midasim/throat_seed1"
)

result$metrics
head(result$fit$results)
```

Set `template` to `"throat"`, `"ibd"`, or `"vaginal"`. When `save_prefix` is
provided, the script writes metrics, taxon-level results, simulation truth, and
the complete result object to files with that prefix.

## Kostic CRC tutorial

An R Jupyter notebook using the Kostic colorectal cancer count data from
`microbiomeMarker` is available at
[`scripts/kostic_data_tutorial.ipynb`](scripts/kostic_data_tutorial.ipynb).
Run it from the repository root so that `devtools::load_all(".")` resolves the
local package. The notebook requires an R kernel and the `microbiomeMarker`,
`phyloseq`, `ggplot2`, `ggrepel`, and `knitr` packages.


## Reference

- Kostic, A. D. et al. (2012). Genomic analysis identifies association of *Fusobacterium* with colorectal carcinoma. *Genome Research*, 22(2), 292–298.