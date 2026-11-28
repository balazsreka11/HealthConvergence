# Health Convergence

Computational Statistics group project analysing mortality convergence between
Hungary and a comparator country using Eurostat cause-of-death data.

## Overview

The project compares age-standardized death rates in 2011 and 2019. A
Kitagawa-style decomposition separates the change in the mortality ratio into:

- Cause-specific rate effects (`Kp`)
- Cause-composition effects (`Kpp`)
- A residual for differences not captured by the selected causes

Parametric bootstrap confidence intervals and waterfall visualizations are
included.

## Files

- `health_group_bootstrap.R`: bootstrap analysis and visualizations.
- `app.R`: interactive Shiny application.


## Run the Analysis

Run `health_group_bootstrap.R` in RStudio to calculate the decomposition and
display the waterfall and confidence-interval plots.

To launch the interactive application:

```r
shiny::runApp("app.R")
```

The first application run downloads the Eurostat data and stores it in
`health_cache.rds`.

## Data and Method

The analysis uses Eurostat table `hlth_cd_asdr2`, with 13 broad ICD-10 causes
and an all-cause category. The bootstrap currently treats standardized rates as
Poisson observations, so its uncertainty estimates should be interpreted with
the limitations. 
