#' Fast 'n' Fair Simultaneous Confidence/Prediction Band for Concurrent Functional Regression
#'
#' @param concurrent_list object returned from running concurrent regression with `fRegress_concurrent`
#' @param interval return a "confidence" or "prediction" band?
#' @param new_dat structure should be a vector if the predictor variables are scalar, otherwise it should be a matrix if functional
#' \itemize{
#'      \item vector length should be Kx1, where K is the number of predictor variables, including the intercept.
#'      \item matrix should be TxK, where T is the number of sampling points and K is the number of predictor variables.
#' }
#' @param time_vec sampling points of the functional observations; default (NULL) sets the sampling points as 1:T, where T is the number of rows in the y matrix returned in the concurrent list
#' @param mse_scalar should the variance estimator be unbiased ("ub") or maximum likelihood ("mle")?
#' @param err_str do the functional errors follow a Gaussian process ("normal") or Student's t process ("t")? "t" is default
#' @param nu0_hat either "singh_df" (default) or a numerical value specifying the degrees of freedom of the Student's t error process
#' @param ... additional arguments that can be passed to `confidence_band`, which creates the final simultaneous band; of note, argument `conf.level` can be set to the desired confidence level for the simultaneous band
#' @return a list containing fast n fair simultaneous confidence or prediction band for the marginal mean of a concurrent functional regression ("simultaneous_band") and an estimate of the degrees of freedom for the Student's t error process ("nu0_hat_est"), if applicable
#' @export
predict_concurrent <- function(concurrent_list, interval = "prediction", new_dat = NULL, time_vec = NULL, mse_scalar = "ub", err_str = "t", nu0_hat = "singh_df", ...) {
  # concurrent_list (list): object returned from running concurrent
  #    regression
  # interval: (char) return confidence or prediction band?
  # new_dat: (vec if scalar covariates, matrix if functional)
  #    x values at which the band should be evaluated;
  #    - if a vector, it should be of length (px1)
  #    - if a matrix, it should be (qxp)
  # time_vec: (vec) sampling points of each curve
  # mse_scalar: (char) which scalar should be used to estimate variance?
  #   -- "mle": 1/n
  #   -- "ub": 1/(n-K)
  #   -- "mm": (n+2)/(n-K)
  # err_str: (char) is the error structure gaussian ("normal"),
  #    or t-distributed ("t")
  # ...: arguments to be passed to confidence_band
  #   --> Could be type, one.sided, upper, grid.size, conf.level,
  #        Bs.sim.size, n_int
  # nu0_hat: (char/numeric) how should the degrees of freedom be estimated?
  #   -- if numeric, that number is used directly as the d.f.
  #   -- if character:
  #       - "singh_df": equations 2.11 and 2.12 from Singh (1988)
  #       - "sut_and_ali": equation 3.10 from Sutradhar and Ali (1986)


  # Check the error structure assumption
  if (!(err_str %in% c("normal", "t"))) {
    stop(
      "Argument err_str (error structure) must be specified",
      "as 'normal' or 't' \n"
    )
  }

  # Get the interval type and number of intervals
  conf_type <- list(...)$type
  # n_ints <- list(...)$n_int

  # Set argument type to correct one, based on 'err_str' argument
  if (is.null(conf_type) & err_str == "normal") {
    conf_type <- "FFSCB.z"
  } else {
    conf_type <- "FFSCB.t"
  }

  # Error structure must also match confidence/prediction band type
  if (err_str == "normal" & !(conf_type == "FFSCB.z")) {
    stop(
      "Argument `err_str` must match match argument `type`. E.g.",
      "if error structure is 'normal', then type must be",
      "'FFSCB.z'."
    )
  } else if (err_str == "t" & !(conf_type == "FFSCB.t")) {
    stop(
      "Argument `err_str` must match match argument `type`. E.g.",
      "if error structure is 't', then type must be",
      "'FFSCB.t'."
    )
  }

  # Extract pieces from model object
  y_hat_mat <- concurrent_list$y_hat_mat
  beta_hat_mat <- concurrent_list$beta_hat_mat
  y_resid_mat <- concurrent_list$y_resid_mat
  y_mat <- concurrent_list$y_mat
  x_array <- concurrent_list$x_array

  # Determine number of observations, number of sampling points,
  # and number of predictor variables
  n_obs <- ncol(y_mat)
  n_sp <- nrow(y_mat)
  n_pred <- dim(x_array)[3]
  mod_df <- n_obs - n_pred

  # Set the time vector
  if (is.null(time_vec)) {
    time_vec <- 1:n_sp
  } else if (length(time_vec) != n_sp) {
    stop(
      "Argument `time_vec` is of different length than",
      " the number of sampling points from the model."
    )
  }

  # Double check that `new_dat` matches the structure of the model
  if (is.vector(new_dat)) {
    if (length(new_dat) != n_pred - 1) {
      stop(
        "Argument `new_dat` must provide an observed value for ",
        "each predictor variable. `new_dat` has a value for ",
        length(new_dat), " predictors, but the model has ",
        n_pred, " predictor variables."
      )
    }
  } else if (is.matrix(new_dat)) {
    if (ncol(new_dat) != n_pred - 1) {
      stop(
        "Argument `new_dat` must provide an observed value for ",
        "each predictor variable. `new_dat` has a value for ",
        length(new_dat), " predictors, but the model has ",
        n_pred, " predictor variables."
      )
    }
  } else {
    stop("Argument `new_dat` must be either a vector or matrix.")
  }


  ## Start constructing the band
  # First: obtain the residuals and estimate MSE(t)
  e_mat <- y_resid_mat

  ### September 6th, 2023: testing out different estimators for the d.f. and
  ### for the variance/covariance
  ## Step 1: d.f. estimation
  nu0_hat_est <- NA
  # if (is.character(nu0_hat)) {
  if (nu0_hat == "singh_df") {
    # Obtain an estimate of quantity 'a'
    a_hat <- ((1 / n_obs) * rowSums(e_mat^4)) /
      (((1 / n_obs) * rowSums(e_mat^2))^2)

    # Drop any estimates less than 1.5, as this will result in negative
    # degrees of freedom, which violates the assumption of df > 4
    a_hat[a_hat < 3] <- NA

    # Using the estimate of a, estimate 'nu_0': df of errors
    nu0_hat_est <- (2 * (2 * a_hat - 3)) / (a_hat - 3)

    # Remove estimates not greater than 4, per assumption, by
    # imputing the values in a sequence with the average of the last
    # non-negative estimate before the sequence, and the first non-
    # negative value after the sequence
    l4_ind <- which(nu0_hat_est <= 4)

    # Remove invalid estimates and infinite values
    nu0_hat_est[l4_ind] <- NA
    nu0_hat_est[is.infinite(nu0_hat_est)] <- NA
    #

    # Take the minimum value of nu0_hat, and replace with a conservative d.f.
    # if none of the estimates were appropriate
    nu0_hat_est <- min(nu0_hat_est, na.rm = T)
    if (is.infinite(nu0_hat_est)) {
      nu0_hat_est <- 4 + .Machine$double.eps
    }
  } else if (is.numeric(nu0_hat)) {
    nu0_hat_est <- nu0_hat
  } else {
    stop("Argument `nu0_hat` must either be numeric or 'singh_df'")
  }

  # } else if (nu0_hat == "sut_and_ali") {
  #   # Estimate nu0_hat
  #   nu0_hat_est <- 2*(3*rowSums(e_mat^2) - (2/n_obs)*sum(colSums(e_mat^4)))/
  #     (3*rowSums(e_mat^2) - (1/n_obs)*sum(colSums(e_mat^4)))
  #
  #   # Remove estimates not greater than 4, per assumption, by
  #   # imputing the values in a sequence with the average of the last
  #   # non-negative estimate before the sequence, and the first non-
  #   # negative value after the sequence
  #   l4_ind <- which(nu0_hat_est <= 4)
  #
  #   # Remove invalid estimates and infinite values
  #   nu0_hat_est[l4_ind] <- NA
  #   nu0_hat_est[is.infinite(nu0_hat_est)] <- NA
  #   #
  #
  #   # Take the minimum value of nu0_hat, and replace with a conservative d.f.
  #   # if none of the estimates were appropriate
  #   nu0_hat_est <- min(nu0_hat_est, na.rm = T)
  #   if (is.infinite(nu0_hat_est)) {
  #     nu0_hat_est <- 4 + .Machine$double.eps
  #   }
  # } else {
  #   stop("Argument `nu0_hat` must either be numeric or one of two character",
  #        "strings: 'singh_df' or 'sut_and_ali'.")
  # }

  ## Step 2: which scalar to use for MSE for estimating variance/covariance?
  if (mse_scalar == "mle") {
    mse_scalar_est <- 1 / n_obs
  } else if (mse_scalar == "ub") {
    mse_scalar_est <- 1 / (n_obs - n_pred)
  } else {
    stop("Arguemtnt `mse_scalar` must be either 'mle' or 'ub'.")
  }
  # } else if (mse_scalar == "mm") {
  #   mse_scalar_est <- (nu0_hat_est - 4)/((nu0_hat_est - 2)*(n_obs - n_pred + 2))
  # }

  ## Step 3: estimate the variance with the things plugged in
  mse_est <- mse_scalar_est * (rowSums(e_mat^2))


  # Second: obtain inverse of (t(X) %*% X)
  x_mat_inv <- array(
    apply(x_array, 1, FUN = function(x) {
      solve(t(x) %*% x)
    }),
    dim = c(n_pred, n_pred, n_sp)
  )


  # Third: using the x vector provided, estimate the mean and cov
  if (is.vector(new_dat)) {
    x_h <- c(1, new_dat)
    cov_x_pc <- unlist(lapply(1:n_sp, FUN = function(i) {
      c((t(x_h) %*% x_mat_inv[, , i] %*% x_h))
    }))

    # Estimate the mean of the predicted value
    mean_est <- c(t(x_h) %*% t(beta_hat_mat))
  } else if (is.matrix(new_dat)) {
    if (nrow(new_dat) != n_sp) {
      stop("Argument `new_dat` must have the same number of rows as
           sampling points.")
    } else {
      X_h <- cbind(rep(1, n_sp), new_dat)
      cov_x_pc <- unlist(lapply(1:n_sp, FUN = function(i) {
        c((t(X_h[i, ]) %*% x_mat_inv[, , i] %*% X_h[i, ]))
      }))

      # Estimate the mean of the predicted value
      mean_est <- c()
      for (i in 1:n_sp) {
        mean_est[i] <- t(X_h[i, ]) %*% beta_hat_mat[i, ]
      }
    }
  } else {
    stop("Argument `new_dat` must be either a vector or matrix")
  }

  if (interval == "confidence") {
    cov_est <- cov_x_pc * mse_est
    # diag(cov_x_pc * mse_est)
  } else if (interval == "prediction") {
    cov_est <- (cov_x_pc + ((nu0_hat_est - 2) / nu0_hat_est)) * mse_est
    # diag((cov_x_pc + ((nu0_hat_est - 2) / nu0_hat_est)) * mse_est)
  }

  # Fourth: estimate the tau function
  tau_est <- ffscb::tau_fun(x = e_mat)

  # Fifth: pass into main confidence band function
  # Note: "confidence" is fixed for int.type because of how that
  #  argument is handled internally with confidence_band and
  # .make_band_FFSCB_t
  simultaneous_band <- confidence_band(
    x = mean_est,
    cov.x = cov_est,
    tau = tau_est,
    df = nu0_hat_est,
    n.curves = n_obs,
    int.type = "confidence",
    type = conf_type,
    ...
  )

  # Finally: return the band and df estimate
  return(list(simultaneous_band, nu0_hat_est))
}
#####################################################################




#' Fast n Fair Simultaneous Confidence Band for the functional parameters of the concurrent regression model
#'
#' @param concurrent_list object returned from running concurrent regression with `fRegress_concurrent`
#' @param interval return a "confidence" or "prediction" band? "confidence" by default
#' @param mse_scalar should the variance estimator be unbiased ("ub") or maximum likelihood ("mle")?
#' @param time_vec sampling points of the functional observations; default (NULL) sets the sampling points as 1:T, where T is the number of rows in the y matrix returned in the concurrent list
#' @param err_str do the functional errors follow a Gaussian process ("normal") or Student's t process ("t")? "t" is default
#' @param nu0_hat either "singh_df" (default) or a numerical value specifying the degrees of freedom of the Student's t error process
#' @param ... additional arguments that can be passed to `confidence_band`, which creates the final simultaneous band; of note, argument `conf.level` (default is 0.95) can be set to the desired confidence level for the simultaneous band
#' @return object `simultaneous_band_list`, which is a list containing a simultaneous confidence/prediction band with `conf.int` confidence for each concurrent regression parameter
#' @export
confint_concurrent <- function(concurrent_list, interval = "confidence", mse_scalar = "ub", time_vec = NULL, err_str = "t", nu0_hat = "singh_df", ...) {
  # concurrent_list (list): object returned from running concurrent
  #    regression
  # interval: (char) return confidence or prediction band?
  # new_dat: (vec if scalar covariates, matrix if functional)
  #    x values at which the band should be evaluated;
  #    - if a vector, it should be of length (px1)
  #    - if a matrix, it should be (qxp)
  # time_vec: (vec) sampling points of each curve
  # mse_loc: (bool) should mse be estimated locally (using only
  # the observations we are trying to predict)
  # err_str: (char) is the error structure gaussian ("normal"),
  #    or t-distributed ("t")
  # ...: arguments to be passed to confidence_band
  #   --> Could be type, one.sided, upper, grid.size, conf.level,
  #        Bs.sim.size, n_int


  # Check the error structure assumption
  if (!(err_str %in% c("normal", "t"))) {
    stop(
      "Argument err_str (error structure) must be specified",
      "as 'normal' or 't' \n"
    )
  }

  # Get the interval type
  conf_type <- list(...)$type

  # Set argument type to correct one, based on 'err_str' argument
  if (is.null(conf_type) & err_str == "normal") {
    conf_type <- "FFSCB.z"
  } else {
    conf_type <- "FFSCB.t"
  }

  # Error structure must also match confidence/prediction band type
  if (err_str == "normal" & !(conf_type %in% c("FFSCB.z"))) {
    stop(
      "Argument `err_str` must match match argument `type`. E.g.",
      "if error structure is 'normal', then type must be",
      "'FFSCB.z'."
    )
  } else if (err_str == "t" & !(conf_type %in% c("FFSCB.t"))) {
    stop(
      "Argument `err_str` must match match argument `type`. E.g.",
      "if error structure is 't', then type must be",
      "'FFSCB.t'."
    )
  }

  # Extract pieces from model object
  y_hat_mat <- concurrent_list$y_hat_mat
  beta_hat_mat <- concurrent_list$beta_hat_mat
  y_resid_mat <- concurrent_list$y_resid_mat
  y_mat <- concurrent_list$y_mat
  x_array <- concurrent_list$x_array

  # Determine number of observations, number of sampling points,
  # and number of predictor variables
  n_obs <- ncol(y_mat)
  n_sp <- nrow(y_mat)
  n_pred <- dim(x_array)[3]
  mod_df <- n_obs - n_pred

  # Set the time vector
  if (is.null(time_vec)) {
    time_vec <- 1:n_sp
  } else if (length(time_vec) != n_sp) {
    stop(
      "Argument `time_vec` is of different length than",
      " the number of sampling points from the model."
    )
  }

  ## Start constructing the band
  # First: obtain the residuals and estimate MSE(t)
  e_mat <- as.matrix(y_mat - y_hat_mat)

  ### September 6th, 2023: testing out different estimators for the d.f. and
  ### for the variance/covariance
  ## Step 1: d.f. estimation
  nu0_hat_est <- NA
  # if (is.character(nu0_hat)) {
  if (nu0_hat == "singh_df") {
    # Obtain an estimate of quantity 'a'
    a_hat <- ((1 / n_obs) * rowSums(e_mat^4)) /
      (((1 / n_obs) * rowSums(e_mat^2))^2)

    # Drop any estimates less than 1.5, as this will result in negative
    # degrees of freedom, which violates the assumption of df > 4
    a_hat[a_hat < 3] <- NA

    # Using the estimate of a, estimate 'nu_0': df of errors
    nu0_hat_est <- (2 * (2 * a_hat - 3)) / (a_hat - 3)

    # Remove estimates not greater than 4, per assumption, by
    # imputing the values in a sequence with the average of the last
    # non-negative estimate before the sequence, and the first non-
    # negative value after the sequence
    l4_ind <- which(nu0_hat_est <= 4)

    # Remove invalid estimates and infinite values
    nu0_hat_est[l4_ind] <- NA
    nu0_hat_est[is.infinite(nu0_hat_est)] <- NA
    #

    # Take the minimum value of nu0_hat, and replace with a conservative d.f.
    # if none of the estimates were appropriate
    nu0_hat_est <- min(nu0_hat_est, na.rm = T)
    if (is.infinite(nu0_hat_est)) {
      nu0_hat_est <- 4 + .Machine$double.eps
    }
  } else if (is.numeric(nu0_hat)) {
    nu0_hat_est <- nu0_hat
  } else {
    stop("Argument `nu0_hat` must either be numeric or 'singh_df'")
  }


  ## Step 2: which scalar to use for MSE for estimating variance/covariance?
  if (mse_scalar == "mle") {
    mse_scalar_est <- 1 / n_obs
  } else if (mse_scalar == "ub") {
    mse_scalar_est <- 1 / (n_obs - n_pred)
  } else {
    stop("Arguemtnt `mse_scalar` must be either 'mle' or 'ub'.")
  }
  # } else if (mse_scalar == "mm") {
  #   mse_scalar_est <- (nu0_hat_est - 4)/((nu0_hat_est - 2)*(n_obs - n_pred + 2))
  # }

  ## Step 3: estimate the variance with the things plugged in
  mse_est <- mse_scalar_est * (rowSums(e_mat^2))

  # Second: obtain inverse of (t(X) %*% X)
  x_mat_inv <- array(
    apply(x_array, 1, FUN = function(x) {
      solve(t(x) %*% x)
    }),
    dim = c(n_pred, n_pred, n_sp)
  )

  c_beta <- matrix(, nrow = n_sp, ncol = n_pred)
  for (i in 1:dim(x_mat_inv)[3]) {
    c_beta[i, ] <- diag(x_mat_inv[, , i])
  }

  # Third: estimate the beta covariance
  cov_list <- vector("list", n_pred)
  for (i in 1:n_pred) {
    if (interval == "confidence") {
      cov_list[[i]] <- diag(c_beta[, i] * mse_est)
    } else if (interval == "prediction") {
      cov_list[[i]] <- diag((c_beta[, i] + 1) * mse_est)
    }
  }

  # Fourth: estimate the tau function
  tau_est <- ffscb::tau_fun(x = e_mat)

  # Fifth: pass into main confidence band function
  # Note: "confidence" is fixed for int.type because of how that
  #  argument is handled internally with confidence_band and
  # .make_band_FFSCB_t
  simultaneous_band_list <- vector("list", n_pred)
  for (i in 1:n_pred) {
    simultaneous_band_list[[i]] <- confidence_band(
      x = beta_hat_mat[, i],
      cov.x = cov_list[[i]],
      tau = tau_est,
      df = nu0_hat,
      n.curves = n_obs,
      int.type = "confidence",
      type = conf_type,
      ...
    )
  }

  # Finally: return the band and df estimate
  return(simultaneous_band_list)
}
#####################################################################
