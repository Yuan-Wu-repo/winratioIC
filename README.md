# winratioIC

`winratioIC` is companion R software for estimating a win ratio with two
prioritized time-to-event outcomes under interval and mixed censoring.

## Supported observation schemes

Version 1.0.0 supports:

1. `margin1 = interval`, `margin2 = interval`
2. `margin1 = right`, `margin2 = interval`

`margin1` is always the higher-priority outcome. The reverse
`interval/right` configuration is not implemented in this version.

## Lower integration bound `tau_l`

The public function includes a user-selected lower integration bound `tau_l`
with

```text
0 <= tau_l < tau
```

- For **interval/interval**, both the primary and secondary win/loss
  contributions are integrated over `[tau_l, tau]`.
- For **right/interval**, the exactly observed/right-censored primary
  contribution remains integrated over `[0, tau]`, while the
  interval-censored secondary contribution is integrated over `[tau_l, tau]`.
- `tau_l = 0` computes the untrimmed exploratory target.
- A strictly positive `tau_l` computes the trimmed target used by the
  positive-lower-bound asymptotic theory, provided the assumptions of the
  methodology are satisfied.

The package does not choose a study-specific positive lower bound automatically.
The analyst should select `tau_l` according to the monitoring-time support and
the estimand/theory used in the analysis.

## Installation

Build and check the package with:

```sh
R CMD build winratioIC
R CMD check winratioIC_1.0.0.tar.gz
```

Install the built source package in R with:

```r
install.packages("winratioIC_1.0.0.tar.gz", repos = NULL, type = "source")
library(winratioIC)
```

## Interval-censored input

For an interval-censored margin, supply lower and upper bounds in `(L,R]`
notation:

- left censored: `L = 0`, finite `R`
- interval censored: `0 < L < R < Inf`
- right censored: finite `L`, `R = Inf`

Example:

```r
margin1 = list(type = "interval", left = "L1", right = "R1")
```

The two interval-censored margins may arise from the same or different visit
processes. No `same_visits` option is required.

Optional `visit_lower` and `visit_upper` column names may be supplied when the
finite inspection pair used for spline knot placement is available. They are
not required for ordinary use.

## Exactly observed / right-censored input

Use `time` and `status`, with `status = 1` for an observed event and
`status = 0` for right censoring:

```r
margin1 = list(type = "right", time = "time1", status = "status1")
```

For the terminal-primary right/interval setting, the secondary interval must
already reflect the information actually observed after any truncation at the
end of primary follow-up.

## Example: two interval-censored margins

```r
library(winratioIC)

dat <- read.csv(system.file(
  "extdata", "example_interval_interval.csv",
  package = "winratioIC"
))

fit <- winratio_ic(
  data = dat,
  treatment = "trt",
  margin1 = list(type = "interval", left = "L1", right = "R1"),
  margin2 = list(type = "interval", left = "L2", right = "R2"),
  tau = 5,
  tau_l = 0.05,
  B = 20,
  k = 6,
  seed = 101
)
fit
```

For the exploratory untrimmed target, use `tau_l = 0`.

## Example: right-censored primary and interval-censored secondary

```r
library(winratioIC)

dat <- read.csv(system.file(
  "extdata", "example_right_interval.csv",
  package = "winratioIC"
))

fit <- winratio_ic(
  data = dat,
  treatment = "trt",
  margin1 = list(type = "right", time = "time1", status = "status1"),
  margin2 = list(type = "interval", left = "L2", right = "R2"),
  tau = 5,
  tau_l = 0.05,
  B = 20,
  k = 6,
  seed = 102
)
fit
```

## Main output

Printing a fitted object reports the censoring scheme, group sizes, analysis
horizon, lower integration bound, whether the target is trimmed or untrimmed,
win and loss probabilities, win ratio, log win ratio, bootstrap standard
error, confidence interval, and bootstrap success count.

The complete fitted object also contains `tau_l`, the estimator label,
per-bootstrap diagnostics, and the two arm-specific fitted spline models.
