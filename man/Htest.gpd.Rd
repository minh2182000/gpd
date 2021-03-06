% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Htest.gpd.R
\name{Htest.gpd}
\alias{Htest.gpd}
\title{Hypothesis test whether a censored data was generated from GPD distribution}
\usage{
Htest.gpd(time, censor, k.hat = NULL, a.hat = NULL)
}
\arguments{
\item{time}{Time to failure or termination}

\item{censor}{Observation/censor indicator (1=observed, 0=right-censored), default to be all observations}

\item{k.hat}{optional, estimate of k, estimated by mle.gpd if either k.hat or a.hat is not given}

\item{a.hat}{optional, estimate of alpha, estimated by mle.gpd if either k.hat or a.hat is not given}
}
\value{
returns a list of the following elements:

\code{Null.Hypo}   States the Null Hypothesis

\code{test.stat}   Test statistic (Y-squared)

\code{df}          Degrees of freedom

\code{p.value}     p-value of the test
}
\description{
\code{Htest.gpd} tests the Null Hypothesis that a right-censored data set follows GPD with the given parameters. Paremeters are estimated with mle.gpd if not given
}
\details{
Methodology was described in:
Pham, M. H., Tsokos, C. P., Choi, B. (2018). Maximum Likelihood Estimation for the Generalized Pareto Distribution and Goodness-Of-Fit Test with Censored Data. Journal of Modern Applied Statistical Methods, 17(2).
}
\examples{
data(melanoma, package = "boot")
status = ifelse(melanoma$status == 1, 1, 0)
Htest.gpd(melanoma$time, status)


}
