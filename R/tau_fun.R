#' This function computes the estimate of the roughness parameter function tau(t) using the pointwise standard deviation of the standardized and differentiated sample functions.
#'
#' @param x Matrix of sample functions (nrow=p, ncol=n, p=number of discretization point, n=sample size).
#' @return tau_t Pointwise standard deviation of the standardized and differentiated sample functions.
#' @examples
#' p <- 200
#' N <- 10
#' rangeval <- c(0, 1)
#' grid <- make_grid(p, rangevals = rangeval)
#' mu <- meanf_poly(grid, params = c(0, 0))
#'
#' # Generate random functions using a stationary
#' # covariance function (homogeneous roughness (HR))
#' cov.m <- make_cov_m(
#'   cov.f = covf_st_matern, grid = grid,
#'   cov.f.params = c(2, 2)
#' )
#' X_HR <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
#'
#' # Generate random functions using non-stationary
#' # covariance function (increasing roughness (IR))
#' cov.m <- make_cov_m(
#'   cov.f = covf_nonst_matern, grid = grid,
#'   cov.f.params = c(3 / 2, 1 / 2, 2)
#' )
#' X_IR <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
#'
#' # Estimate tau(t):
#' tau_HR <- tau_fun(X_HR)
#' tau_IR <- tau_fun(X_IR)
#'
#' # Plot data and estimated tau() functions
#' par(mfrow = c(2, 2))
#' matplot(
#'   x = grid, y = X_HR, type = "l", main = "Homogeneous Roughness",
#'   ylab = "X(t)", xlab = ""
#' )
#' matplot(
#'   x = grid, y = X_IR, type = "l", main = "Increasing Roughness",
#'   ylab = "X(t)", xlab = ""
#' )
#' plot(
#'   x = grid, y = tau_HR, type = "l", main = "Homogeneous Roughness",
#'   ylab = "tau(t)", xlab = "", ylim = range(tau_HR, tau_IR)
#' )
#' plot(
#'   x = grid, y = tau_IR, type = "l", main = "Increasing Roughness",
#'   ylab = "tau(t)", xlab = "", ylim = range(tau_HR, tau_IR)
#' )
#' par(mfrow = c(1, 1))
#' @export tau_fun
tau_fun <- function(x) {
  ##
  x_scl <- t(apply(x, 1, scale)) # pointwise standardization (mean=0, sd=1)
  x_scl <- x_scl
  x_p <- apply(x_scl, 2,
    FUN = function(yy) { # differentiation
      xx <- seq(0, 1, len = length(yy))
      fn <- stats::splinefun(x = xx, y = yy, method = "natural")
      pracma::fderiv(f = fn, x = xx, n = 1, h = diff(xx)[1], method = "central")
    }
  )
  tau_t <- apply(x_p, 1, stats::sd) # pointwise sd
  return(tau_t)
}




#' R-function for computing tau from fragmentary functional data.
#' Caution: only one single fragment per function is assumed.
#'
#' @param X_mat Matrix of sample functions (nrow=p, ncol=n, p=number of discretization point, n=sample size).
#' @param grid_mat Matrix (pxn) of grid points over which the fragments are observed.
# @export
tau_fragments <- function(X_mat, grid_mat) {
  ##
  ## Caution: this function assumes *single* fragments per function
  ##
  ## Assumed dimensions:
  ## ncol == sample size (i.e., number of functions)
  ## nrow == number of discretization points
  ##
  ## Check [0,1] restriction
  if (!all(range(grid_mat, na.rm = TRUE) == c(0, 1))) {
    stop("grid_mat has to be within [0,1]")
  }
  ##
  X_scl <- t(apply(X_mat, 1, scale))
  ## Derivative
  X_p <- matrix(NA, nrow = nrow(X_mat), ncol = ncol(X_mat))
  for (j in 1:ncol(X_p)) {
    fn <- stats::splinefun(x = grid_mat[, j], y = X_scl[, j], method = "natural")
    X_p[!is.na(grid_mat[, j]), j] <- pracma::fderiv(
      f = fn, x = c(stats::na.omit(grid_mat[, j])), n = 1,
      h = diff(c(stats::na.omit(grid_mat[, j])))[1], method = "central"
    )
  }
  ##
  tau_t <- apply(X_p, 1, function(x) stats::sd(x, na.rm = TRUE))
  return(tau_t)
}





#' This function computes the estimate of the roughness parameter function tau(t) using the covariance function (given as a matrix) of the functional data.
#'
#' @param cov_mat Matrix (pxp) of evaluated covariance function (p=number of discretization point). Caution: It is assumed that the evaluation grid is within [0,1].
#' @param warn Option for printing warnings
#' @return tau_t Estimate of the roughness parameter function tau(t)
#' @examples
#' p <- 200
#' N <- 10
#' rangeval <- c(0, 1)
#' grid <- make_grid(p, rangevals = rangeval)
#' mu <- meanf_poly(grid, params = c(0, 0))
#'
#' # Generate random functions using a stationary
#' # covariance function (homogeneous roughness (HR))
#' cov.m <- make_cov_m(
#'   cov.f = covf_st_matern, grid = grid,
#'   cov.f.params = c(2, 2, 2)
#' )
#' X_HR <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
#'
#' # Generate random functions using non-stationary
#' # covariance function (increasing roughness (IR))
#' cov.m <- make_cov_m(
#'   cov.f = covf_nonst_matern, grid = grid,
#'   cov.f.params = c(3 / 2, 1 / 2, 2)
#' )
#' X_IR <- make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
#'
#' # Estimate covariance functions
#' hat_mu_HR <- rowMeans(X_HR)
#' hat_cov_HR <- crossprod(t(X_HR - hat_mu_HR)) / (N - 1)
#'
#' hat_mu_IR <- rowMeans(X_IR)
#' hat_cov_IR <- crossprod(t(X_IR - hat_mu_IR)) / (N - 1)
#'
#' # Estimate tau(t):
#' tau_HR <- cov2tau_fun(hat_cov_HR)
#' tau_IR <- cov2tau_fun(hat_cov_IR)
#'
#' # Plot data and estimated tau() functions
#' par(mfrow = c(2, 2))
#' matplot(
#'   x = grid, y = X_HR, type = "l", main = "Homogeneous Roughness",
#'   ylab = "X(t)", xlab = ""
#' )
#' matplot(
#'   x = grid, y = X_IR, type = "l", main = "Increasing Roughness",
#'   ylab = "X(t)", xlab = ""
#' )
#' plot(
#'   x = grid, y = tau_HR, type = "l", main = "Homogeneous Roughness",
#'   ylab = "tau(t)", xlab = "", ylim = range(tau_HR, tau_IR)
#' )
#' plot(
#'   x = grid, y = tau_IR, type = "l", main = "Increasing Roughness",
#'   ylab = "tau(t)", xlab = "", ylim = range(tau_HR, tau_IR)
#' )
#' par(mfrow = c(1, 1))
#' @export cov2tau_fun
cov2tau_fun <- function(cov_mat, warn = FALSE) {
  ## 'cov_mat' denotes the sample covariance function computed from the sample functions X_1(t),...X_n(t)
  ## Caution: assumed grid is in [0,1]
  ##
  p <- ncol(cov_mat)
  grid <- seq(0, 1, len = p)
  corr_mat <- stats::cov2cor(cov_mat)
  ## computing the numeric approximation to c_12(t,s) := \partial^2 c(t,s)/(\partial t \partial t)
  ## with c_12(t,s) evaluated at t=t and s=t
  a1 <- 1 # Equals 1: corr_mat[cbind(2:p,2:p)]        # corr(t+h, t+h) with h=0.005, and t\in{0.005,0.015,0.025,...,0.995}
  a2 <- corr_mat[cbind(1:(p - 1), 2:p)] # corr(t-h, t+h)
  a3 <- 1 # Equals 1: corr_mat[cbind(1:(p-1),1:(p-1))]# corr(t-h, t-h)
  ##
  h <- diff(grid)[1] / 2
  suppressWarnings(
    tau <- sqrt(c(a1 - 2 * a2 + a3) / (4 * h^2))
  )
  if (any(is.na(tau))) {
    if (warn) {
      warning("\nCaution there were negative tau-values!\n
             Probably, this was caused by a cov-estimate\n
             computed from fragmentary functions.\n
             We use only the positive values and fill the \n
             missing tau-values by approximations using linear interpolations.")
    }
    if (is.na(tau[1])) {
      tau[1] <- utils::head(c(stats::na.omit(tau)), n = 1)
    }
    if (is.na(tau[p - 1])) {
      tau[p - 1] <- utils::tail(c(stats::na.omit(tau)), n = 1)
    }
    tau <- stats::approx(x = c(1:(p - 1))[!is.na(tau)], y = tau[!is.na(tau)], xout = c(1:(p - 1)))$y
  }
  ## tau has length p-1, so we interpolate to get length p
  xx <- c(0, grid[-p] + h, 1)
  yy <- c(tau[1], tau, tau[length(tau)])
  tau_t <- stats::approx(y = yy, x = xx, xout = grid)$y
  ##
  return(tau_t)
}

# # This function computes the estimate of the roughness parameter function tau(t) using the pointwise standard deviation of the standardized and differentiated sample functions.
# #
# # @param x Matrix of sample functions (nrow=p, ncol=n, p=number of discretization point, n=sample size).
# # @param df Possibility to take into account the degrees of freedom (typically df=N-1). If df is provided, the internally standardized sample functions are scaled by ((df-2)/df)^(1/2). If df=null (default), there will be no scaling.
# # @return tau_t Pointwise standard deviation of the standardized and differentiated sample functions.
# # @examples
# # p         <- 200
# # N         <- 10
# # rangeval  <- c(0,1)
# # grid      <- make_grid(p, rangevals=rangeval)
# # mu        <- meanf_poly(grid, params = c(0,0))
# #
# # # Generate random functions using a stationary
# # # covariance function (homogeneous roughness (HR))
# # cov.m = make_cov_m(cov.f = covf_st_matern, grid=grid,
# # cov.f.params=c(2,2,2))
# # X_HR  <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
# #
# # # Generate random functions using non-stationary
# # # covariance function (increasing roughness (IR))
# # cov.m = make_cov_m(cov.f = covf_st_matern.warp.power, grid=grid,
# # cov.f.params=c(1, 1, 1, 3))
# # X_IR  <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
# #
# # # Estimate tau(t):
# # tau_HR  <- tau_fun(X_HR, df=N-1)
# # tau_HR2 <- tau_fun(X_HR)
# # tau_IR  <- tau_fun(X_IR, df=N-1)
# # tau_IR2 <- tau_fun(X_IR)
# #
# # # Plot data and estimated tau() functions
# # par(mfrow=c(2,2))
# # matplot(x=grid, y=X_HR, type="l", main="Homogeneous Roughness",
# # ylab="X(t)", xlab="")
# # matplot(x=grid, y=X_IR, type="l", main="Increasing Roughness",
# # ylab="X(t)", xlab="")
# # matplot(x=grid, y=cbind(tau_HR,tau_HR2),  type="l", main="Homogeneous Roughness",
# # ylab="tau(t)", xlab="", ylim=range(tau_HR, tau_IR))
# # matplot(x=grid, y=cbind(tau_IR,tau_IR2),  type="l", main="Increasing Roughness",
# # ylab="tau(t)", xlab="", ylim=range(tau_HR, tau_IR))
# # par(mfrow=c(1,1))
# # @export tau_fun
# tau_fun <- function(x, df=NULL){
#   ##
#   if( is.null(df)       ){
#     scl <- 1
#   }else{
#     if(df>2){ scl <- (df-2)/df }else{ scl <- 1 }
#   }
#   ##
#   x_scl    <- t(apply(x, 1, scale)) # pointwise standardization (mean=0, sd=1)
#   x_scl    <- x_scl * sqrt(scl)
#   x_p      <- apply(x_scl, 2,
#                     FUN=function(yy){     # differentiation
#                       xx <- seq(0,1,len=length(yy))
#                       fn <- stats::splinefun(x = xx, y = yy, method = "natural")
#                       pracma::fderiv(f = fn, x = xx, n = 1, h = diff(xx)[1], method = "central")})
#   tau_t    <- apply(x_p, 1, stats::sd)           # pointwise sd
#   return(tau_t)
# }
#
#
#
#
#
# # This function computes the estimate of the roughness parameter function tau(t) using the covariance function (given as a matrix) of the functional data.
# #
# # @param cov_mat Matrix (pxp) of evaluated covariance function (p=number of discretization point). Caution: It is assumed that the evaluation grid is within [0,1].
# # @param df Possibility to take into account the degrees of freedom (typically df=N-1). If df is provided, the internally computed correlation function is scaled by (df-2)/df. If df=null (default), there will be no scaling.
# # @return tau_t Estimate of the roughness parameter function tau(t)
# # @examples
# # p         <- 200
# # N         <- 10
# # rangeval  <- c(0,1)
# # grid      <- make_grid(p, rangevals=rangeval)
# # mu        <- meanf_poly(grid, params = c(0,0))
# #
# # # Generate random functions using a stationary
# # # covariance function (homogeneous roughness (HR))
# # cov.m = make_cov_m(cov.f = covf_st_matern, grid=grid,
# # cov.f.params=c(2,2,2))
# # X_HR  <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
# #
# # # Generate random functions using non-stationary
# # # covariance function (increasing roughness (IR))
# # cov.m = make_cov_m(cov.f = covf_st_matern.warp.power, grid=grid,
# # cov.f.params=c(1, 1, 1, 3))
# # X_IR  <-  make_sample(mean.v = mu, cov.m = cov.m, N = N, dist = "rnorm")
# #
# # # Estimate covariance functions
# # hat_mu_HR   <- rowMeans(X_HR)
# # hat_cov_HR  <- crossprod(t(X_HR - hat_mu_HR)) / (N-1)
# #
# # hat_mu_IR   <- rowMeans(X_IR)
# # hat_cov_IR  <- crossprod(t(X_IR - hat_mu_IR)) / (N-1)
# #
# # # Estimate tau(t):
# # tau_HR  <- cov2tau_fun(hat_cov_HR, df=N-1)
# # tau_HR2 <- cov2tau_fun(hat_cov_HR)
# # tau_IR  <- cov2tau_fun(hat_cov_IR, df=N-1)
# # tau_IR2 <- cov2tau_fun(hat_cov_IR)
# #
# # # Plot data and estimated tau() functions
# # par(mfrow=c(2,2))
# # matplot(x=grid, y=X_HR, type="l", main="Homogeneous Roughness",
# # ylab="X(t)", xlab="")
# # matplot(x=grid, y=X_IR, type="l", main="Increasing Roughness",
# # ylab="X(t)", xlab="")
# # matplot(x=grid, y=cbind(tau_HR,tau_HR2),  type="l", main="Homogeneous Roughness",
# # ylab="tau(t)", xlab="", ylim=range(tau_HR, tau_IR))
# # matplot(x=grid, y=cbind(tau_IR,tau_IR2),  type="l", main="Increasing Roughness",
# # ylab="tau(t)", xlab="", ylim=range(tau_HR, tau_IR))
# # par(mfrow=c(1,1))
# # @export cov2tau_fun
# cov2tau_fun <- function(cov_mat, df=NULL){
#   ## 'cov_mat' denotes the sample covariance function computed from the sample functions X_1(t),...X_n(t)
#   ## Caution: assumed grid is in [0,1]
#   ##
#   if( is.null(df)       ){
#     scl <- 1
#   }else{
#     if(df>2){ scl <- (df-2)/df }else{ scl <- 1 }
#   }
#   ##
#   p        <- ncol(cov_mat)
#   grid     <- seq(0,1,len=p)
#   corr_mat <- stats::cov2cor(cov_mat) * scl
#   ## computing the numeric approximation to c_12(t,s) := \partial^2 c(t,s)/(\partial t \partial t)
#   ## with c_12(t,s) evaluated at t=t and s=t
#   a1  <- corr_mat[cbind(2:p,2:p)]        # corr(t+h, t+h) with h=0.005, and t\in{0.005,0.015,0.025,...,0.995}
#   a2  <- corr_mat[cbind(1:(p-1),2:p)]    # corr(t-h, t+h)
#   a3  <- corr_mat[cbind(1:(p-1),1:(p-1))]# corr(t-h, t-h)
#   ##
#   h   <- diff(grid)[1]/2
#   tau <- sqrt( c(a1- 2*a2 + a3) /( 4* h^2 ) )
#   ## tau has length p-1, so we interpolate to get length p
#   xx    <- c(0,      grid[-p] + h, 1               )
#   yy    <- c(tau[1], tau,          tau[length(tau)])
#   tau_t <- stats::spline(y=yy, x=xx, xout = grid, method = "natural")$y
#   ##
#   return(tau_t)
# }
