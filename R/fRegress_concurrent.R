#' Functional Concurrent Regression without any smoothing or regularization
#'
#' @param y_mat matrix of response curves, with number of rows equal to the number of sampling points, T, and the number of columns equal to the number of observations, N (T x N)
#' @param x_array first two dimensions are the same as `y_mat`, while the third dimension represents the number of predictor variables, K (T x N x K)
#' @param intercept logical argument: does `x_array` contain a matrix of one's representing the intercept? `FALSE` by default, so the function creates the intercept term
#' \itemize{
#'      \item Note: if there's only one predictor variable aside from the intercept, `x_array` can be a T x N matrix with `intercept` set to `FALSE`
#' }
#' @return a list containing the predicted functional observations, `y_hat_mat`, the estimated functional beta parameters, `beta_hat_mat`, the functional residuals, `y_resid_mat`, the original functional observations, `y_mat`, the original functional/scalar predictor variables, `x_array`, estimated functional DFFITS, `dffits_mat`, and estimated functional studentized residuals, `ext_stu_res_mat`
#' @export
fRegress_concurrent <- function(y_mat, x_array, intercept = FALSE) {
  # y_mat (mat): matrix of response curves, with nrows = number of
  #              sampling points and ncols = number of observations
  # x_array (array): first two dimensions are same as y_mat, while
  #                  the third dimension represents number of pred's
  # intercept (bool): did you include the intercept column in X? T/F

  # Check structure of x_array and set number of preds accordingly
  if (is.matrix(x_array)) {
    n_pred <- 1
  } else if (is.array(x_array)) {
    n_pred <- dim(x_array)[3]
  } else {
    stop("Wrong structure for x_array argument.")
  }

  # Determine number of observations, number of sampling points,
  # and number of predictor variables
  n_obs <- ncol(y_mat)
  n_sp <- nrow(y_mat)

  if (intercept) {
    # Need empty data frames for predicted values, estimated beta
    # coefficients, and each residual
    y_hat_mat <- matrix(, nrow = n_sp, ncol = n_obs)
    beta_hat_mat <- matrix(, nrow = n_sp, ncol = n_pred)
    y_resid_mat <- matrix(, nrow = n_sp, ncol = n_obs)

    # Loop over the sampling points, fitting a linear regression at
    # each sampling point
    for (i in 1:n_sp) {
      # Objects for regression
      ys <- y_mat[i, ]

      # Set based on x_array structure
      if (is.matrix(x_array)) {
        x_mat <- x_array[i, ]
        mod_form_temp <- stats::as.formula(paste0("ys ~ -1 + x_mat"))
      } else {
        x_mat <- x_array[i, , ]
        mod_form_temp <- stats::as.formula(paste0(
          "ys ~ -1 + ",
          paste0("x_mat[,", 1:n_pred, "]", collapse = "+")
        ))
      }

      # SLR
      slr_out <- stats::lm(formula = mod_form_temp)

      # Save the values
      y_hat_mat[i, ] <- slr_out$fitted.values
      beta_hat_mat[i, ] <- slr_out$coefficients
      y_resid_mat[i, ] <- slr_out$residuals
    }
  } else {
    # Need empty data frames for predicted values, estimated beta
    # coefficients, and each residual
    y_hat_mat <- matrix(, nrow = n_sp, ncol = n_obs)
    beta_hat_mat <- matrix(, nrow = n_sp, ncol = n_pred + 1)
    # Note: 1 is added above for the intercept term
    y_resid_mat <- matrix(, nrow = n_sp, ncol = n_obs)
    ext_stu_res_mat <- matrix(, nrow = n_sp, ncol = n_obs)
    dffits_mat <- matrix(, nrow = n_sp, ncol = n_obs)

    # Loop over the sampling points, fitting a linear regression at
    # each sampling point
    for (i in 1:n_sp) {
      # Objects for regression
      ys <- y_mat[i, ]

      # Set based on x_array structure
      if (is.matrix(x_array)) {
        x_mat <- x_array[i, ]
        mod_form_temp <- stats::as.formula(paste0("ys ~ x_mat"))
      } else {
        x_mat <- x_array[i, , ]
        mod_form_temp <- stats::as.formula(paste0(
          "ys ~ ", paste0("x_mat[,", 1:n_pred, "]", collapse = "+")
        ))
      }

      # SLR
      slr_out <- stats::lm(formula = mod_form_temp)

      # Save the values
      y_hat_mat[i, ] <- slr_out$fitted.values
      beta_hat_mat[i, ] <- slr_out$coefficients
      y_resid_mat[i, ] <- slr_out$residuals
      dffits_mat[i, ] <- stats::dffits(slr_out)
      ext_stu_res_mat[i, ] <- MASS::studres(slr_out)
    }

    # Add a matrix of ones for the intercept to the x_array
    x_array <- array(
      data = c(matrix(1, nrow = n_sp, ncol = n_obs), x_array),
      dim = c(nrow = n_sp, ncol = n_obs, n_pred + 1)
    )
  }

  # Return back the fitted objects
  return(list(
    y_hat_mat = y_hat_mat,
    beta_hat_mat = beta_hat_mat,
    y_resid_mat = y_resid_mat,
    y_mat = y_mat,
    x_array = x_array,
    dffits_mat = dffits_mat,
    ext_stu_res_mat = ext_stu_res_mat
  ))
}
#####################################################################
