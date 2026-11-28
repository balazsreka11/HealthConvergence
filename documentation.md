# Health Convergence Analysis

## Purpose

This project (as a group projec for Computational Statistics course) studies whether Hungary's age-standardized mortality profile
converged toward Denmark's between 2011 and 2019. It compares cause-specific
death rates and decomposes the change into rate and composition effects.

## Data

The analysis uses Eurostat table `hlth_cd_asdr2`: **standardized death rates by
cause of death**.

The following settings are used:

- Countries: Hungary (`HU`) and Denmark (`DK`)
- Years: 2011 and 2019
- Sex: total population
- Age: total population
- Unit: age-standardized death rate per 100,000

The data contains 13 broad ICD-10 cause groups plus an all-cause category:

| Code | Cause group |
|---|---|
| `A_B` | Infectious and parasitic diseases |
| `C00-D48` | Neoplasms |
| `D50-D89` | Blood and immune disorders |
| `E` | Endocrine and metabolic diseases |
| `F` | Mental and behavioural disorders |
| `G_H` | Nervous system and sense organs |
| `I` | Circulatory diseases |
| `J` | Respiratory diseases |
| `K` | Digestive diseases |
| `L` | Skin diseases |
| `M` | Musculoskeletal diseases |
| `N` | Genitourinary diseases |
| `V01-Y89` | External causes |

Eurostat changed the all-cause code between the two years. `TOTAL` is used for
2011 and `A-R_V-Y` for 2019; both are harmonized as `ALL_CAUSES`.

## Data Structure

After filtering and reshaping, each cause has four observations:

| Variable | Meaning |
|---|---|
| `inv_0` | Hungary, 2011 |
| `inv_1` | Hungary, 2019 |
| `comp_0` | Denmark, 2011 |
| `comp_1` | Denmark, 2019 |

The values are **age-standardized rates**, not raw death counts.

## Convergence Measure

For each cause, the Hungary-to-Denmark mortality ratio is calculated as:

```text
V_0 = Hungary 2011 / Denmark 2011
V_1 = Hungary 2019 / Denmark 2019
k   = V_1 - V_0
```

Interpretation:

- `V = 1`: equal mortality rates.
- `V > 1`: Hungary has a higher rate.
- `V < 1`: Hungary has a lower rate.
- A negative `k` indicates convergence when Hungary initially has the higher
  rate.

## Kitagawa Decomposition

The 13 individual causes are used to split the ratio change into two effects.

### Rate Effect (`Kp`)

`Kp` measures the contribution of a change in the cause-specific
Hungary-to-Denmark ratio.

- `Kp < 0`: contribution toward convergence.
- `Kp > 0`: contribution toward divergence.

`sum_Kp` is the total rate effect.

### Composition Effect (`Kpp`)

`Kpp` measures the effect of changes in the distribution of Denmark's mortality
across causes. It is based on each cause's share of Denmark's summed mortality
rates.

The sum of the cause-specific values is reported as `composition_effect`.

### Residual

The 13 selected cause rates may not reproduce the published all-cause rate
exactly. The remaining difference is therefore displayed explicitly:

```text
Residual = observed all-cause change - sum_Kp - composition_effect
```

## Bootstrap

A parametric bootstrap creates 10,000 simulated observations for every
country-year-cause rate using:

```text
Simulated rate ~ Poisson(observed rate)
```

The 2.5th and 97.5th percentiles provide 95% confidence intervals for:

- Mortality rates
- Ratio changes (`k`)
- Rate effects (`Kp`)
- Composition effects (`Kpp`)
- Total decomposition effects

**Limitation:** the Poisson model is more naturally applied to death counts than
to standardized rates. A more precise analysis would require death counts,
population exposure, or standardization weights.

## Plots

### Waterfall Chart

The waterfall chart shows:

```text
2011 all-cause ratio
-> 13 cause-specific rate contributions (Kp)
-> total composition effect
-> residual
-> 2019 all-cause ratio
```

- Grey bars: observed ratios at the beginning and end.
- Green bars: contributions toward convergence.
- Red bars: contributions toward divergence.
- Gold bar: residual.

This plot identifies which causes contributed most to the overall change.

### Effect and Confidence Interval Plot

This plot displays the cause-specific `Kp` and `Kpp` estimates with their 95%
bootstrap intervals in separate panels.

The dashed vertical line represents zero:

- Intervals below zero indicate convergence contributions.
- Intervals above zero indicate divergence contributions.
- Intervals crossing zero indicate an uncertain direction.

Confidence intervals are shown separately from the waterfall because simulated
cause contributions are correlated and do not behave as independent waterfall
steps.

## Main Outputs

The analysis produces:

- `rates`: cleaned and reshaped mortality rates.
- `convergence`: `V_0`, `V_1`, and `k`.
- `kitagawa_results`: cause-specific `Kp` and `Kpp` values.
- `kitagawa_summary`: total rate and composition effects.
- `ci_rates`: mortality-rate confidence intervals.
- `ci_kitagawa`: decomposition confidence intervals by cause.
- `ci_totals`: confidence intervals for the total effects.
- `waterfall_plot`: decomposition from 2011 to 2019.
- `effect_ci_plot`: cause effects with bootstrap uncertainty.
