# CABNB

Compositionality-aware, bias-corrected negative binomial model for differential abundance analysis


## Installation

```R
devtools::install_github("noriakis/CABNB")
```

## Example usage

```R
data(cabnb_example_counts)
data(cabnb_example_metadata)
fit <- cabnb_fit(
  cabnb_example_counts,
  cabnb_example_metadata,
  formula = ~ group,
  coef = "grouptreatment"
)
plot(fit)
```

## Run one MIDASim condition

Install MIDASim and the  PRROC package, load CABNB, and source
the benchmark function.

```R
source("scripts/midasim_one_condition.R")

bench <- cabnb_run_midasim_one_condition()
bench$metrics
```