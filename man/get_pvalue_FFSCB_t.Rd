% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_pvalue.R
\name{get_pvalue_FFSCB_t}
\alias{get_pvalue_FFSCB_t}
\title{FFSCB p-value (t-distr)}
\usage{
get_pvalue_FFSCB_t(
  x,
  x0 = NULL,
  tau,
  diag.cov,
  df,
  eval.points = NULL,
  n_int = 5
)
}
\arguments{
\item{x}{Functional parameter estimate (for instance, the empirical mean function).}

\item{x0}{Functional parameter under the null hypothesis. Default: zero.}

\item{tau}{Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().}

\item{diag.cov}{The diagonal of Cov(x), in which x is the functional estimator. For instance, the diagonal of the discretized covariance function of the empirical mean function x.}

\item{df}{Degrees of freedom parameter for the t-distribution based band 'FFSCB.t'. (Typically, df=N-1)}

\item{eval.points}{Evaluation points (in [0,1]) at which the pvalues should be computed.}

\item{n_int}{Number of intervals parameter used by the function make_band_FFSCB_t()}
}
\description{
FFSCB p-value (t-distr)
}
\examples{
# Generate a sample
p <- 200
N <- 80
grid <- make_grid(p, rangevals = c(0, 1))
mu0 <- meanf_poly(grid, c(0, 1))
names(mu0) <- grid
mu <- meanf_poly(grid, c(0, 1.1))
names(mu) <- grid
cov.m <- make_cov_m(cov.f = covf_st_matern, grid = grid, cov.f.params = c(2 / 2, 1, 1))
sample <- make_sample(mu, cov.m, N)

# Compute the estimate and its covariance
hat.mu <- rowMeans(sample)
hat.cov <- crossprod(t(sample - hat.mu)) / N
hat.cov.mu <- hat.cov / N
hat.tau <- tau_fun(sample)

# Compute simultaneous pvalue function
pval <- get_pvalue_FFSCB_t(
  x = hat.mu, x0 = mu0, tau = hat.tau,
  diag.cov = diag(hat.cov.mu), df = N - 1,
  eval.points = c(0.25, 0.75)
)
pval
}
\references{
Liebl, D. and Reimherr, M. (2022+). Fast and fair simultaneous confidence bands.
}
