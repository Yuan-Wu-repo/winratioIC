\name{winratio_ic}
\alias{winratio_ic}
\alias{print.winratio_ic}
\title{Win Ratio Estimation with Interval and Mixed Censoring}
\description{
Estimates the win probability, loss probability, win ratio, log win ratio,
and weighted-bootstrap uncertainty for two prioritized time-to-event outcomes.
The higher-priority outcome is always supplied as \code{margin1}.
}
\usage{
winratio_ic(
  data,
  treatment,
  margin1,
  margin2,
  tau,
  tau_l = 0,
  B = 200L,
  k = 6L,
  parallel = FALSE,
  seed = NULL,
  conf_level = 0.95
)

\method{print}{winratio_ic}(x, ...)
}
\arguments{
  \item{data}{A data frame with one row per subject.}
  \item{treatment}{Name of the treatment-group column. Numeric treatment must
  be coded 0 for control and 1 for treatment. A two-level factor is also
  accepted, with the first observed factor level treated as control and the
  second as treatment. Character labels \code{"Control"} and
  \code{"Treatment"} are recognized case-insensitively.}
  \item{margin1}{A list describing the higher-priority outcome. For an
  interval-censored margin use
  \code{list(type="interval", left="L1", right="R1")}. For an exactly
  observed/right-censored margin use
  \code{list(type="right", time="time1", status="status1")}.}
  \item{margin2}{A list describing the secondary outcome. In the current
  version the supported configurations are interval/interval and
  right/interval.}
  \item{tau}{Positive finite analysis horizon.}
  \item{tau_l}{Lower integration bound satisfying \code{0 <= tau_l < tau}.
  For interval/interval, both prioritized win/loss contributions are
  integrated over \eqn{[\tau_l,\tau]}. For right/interval, the primary
  contribution remains over \eqn{[0,\tau]} and the interval-censored
  secondary contribution uses \eqn{[\tau_l,\tau]}. \code{tau_l=0} computes
  the untrimmed exploratory target.}
  \item{B}{Number of weighted-bootstrap multiplier draws. Must be at least 2.}
  \item{k}{Spline-sieve dimension. Default is 6.}
  \item{parallel}{Logical. Parallel execution is not implemented in version
  1.0.0 and must be \code{FALSE}.}
  \item{seed}{Optional random-number seed used for the bootstrap multiplier
  draws.}
  \item{conf_level}{Confidence level for the normal interval on the log win
  ratio scale, transformed back to the win-ratio scale.}
  \item{x}{An object returned by \code{winratio_ic}.}
  \item{...}{Additional arguments (currently unused).}
}
\details{
For an interval-censored margin, observations use the usual \eqn{(L,R]}
representation:
\itemize{
  \item left censored: \code{L = 0} and finite \code{R};
  \item interval censored: \code{0 < L < R < Inf};
  \item right censored: finite \code{L} and \code{R = Inf}.
}
Exact observations with \code{left == right} are not supported under
\code{type="interval"}.

For \code{type="right"}, \code{status=1} denotes an observed event and
\code{status=0} denotes right censoring.

The lower bound \code{tau_l} changes the target win/loss functional, not the
likelihood fit. In the interval/interval setting, both prioritized components
use the restricted interval \eqn{[\tau_l,\tau]}. In the right/interval
setting, the exactly observed/right-censored primary comparison remains on
\eqn{[0,\tau]}, while only the interval-censored secondary contribution is
restricted to \eqn{[\tau_l,\tau]}. The same \code{tau_l} is used for the
point estimate and every weighted-bootstrap replicate.

A strictly positive lower bound is the setting covered by the
positive-lower-bound asymptotic theory, subject to the assumptions of the
methodology. Setting \code{tau_l=0} is allowed to compute the corresponding
untrimmed exploratory estimand but is not covered by that positive-lower-bound
asymptotic result.

When the interval inspection pair used by the original spline-knot rule is
available, an interval margin may additionally specify \code{visit_lower} and
\code{visit_upper}. These optional columns are used only for knot placement;
the censoring category remains determined by \code{left} and \code{right}.
Ordinary users do not need these optional columns.

In the right/interval terminal-outcome setting, the supplied interval for the
secondary outcome must already represent the observed secondary information
after any truncation at the end of primary follow-up. The function does not
reconstruct visit histories or re-apply study-specific truncation rules.

The reverse interval/right configuration is not implemented in this version.
}
\value{
An object of class \code{"winratio_ic"}. Important components include:
\itemize{
  \item \code{tau_l}: lower integration bound used for the functional;
  \item \code{estimator}: \code{"trimmed"} when \code{tau_l>0} and
  \code{"untrimmed"} when \code{tau_l=0};
  \item \code{win_probability}: estimated probability that treatment wins;
  \item \code{loss_probability}: estimated probability that treatment loses;
  \item \code{win_ratio}: estimated win/loss ratio;
  \item \code{log_win_ratio}: logarithm of the win ratio;
  \item \code{se_log_win_ratio}: bootstrap standard error of the log win ratio;
  \item \code{ci_win_ratio}: confidence interval on the win-ratio scale;
  \item \code{bootstrap}: bootstrap success counts and diagnostics;
  \item \code{fit}: fitted control- and treatment-arm spline models.
}
}
\author{
Zhiguo Li \email{zhiguo.li@duke.edu}
}
\examples{
\dontrun{
## Interval / interval: positive lower bound
fit <- winratio_ic(
  data = dat,
  treatment = "trt",
  margin1 = list(type = "interval", left = "L1", right = "R1"),
  margin2 = list(type = "interval", left = "L2", right = "R2"),
  tau = 5,
  tau_l = 0.05,
  B = 200,
  k = 6,
  seed = 123
)
fit

## Right / interval: primary remains [0,tau], secondary uses [tau_l,tau]
fit2 <- winratio_ic(
  data = dat2,
  treatment = "trt",
  margin1 = list(type = "right", time = "time1", status = "status1"),
  margin2 = list(type = "interval", left = "L2", right = "R2"),
  tau = 5,
  tau_l = 0.05,
  B = 200,
  k = 6,
  seed = 123
)
fit2

## Exploratory untrimmed target: set tau_l = 0.
}
}
