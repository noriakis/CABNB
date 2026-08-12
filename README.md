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
