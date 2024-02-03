#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### August 16th, 2023                                             ###
###                                                               ###
###   Updated simulation study. This simulation compares the Fast ###
### n Fair simultaneous confidence/prediction bands for concurrent###
### functional regression by Creutzinger, Liebl, and Sharp (2024+)###
### versus the conformal inference bands of Diquigiovanni et al.  ###
### (2022). There is also functionality to recreate the plots in  ###
### Figures 2 and 3.                                              ###
#####################################################################


## Starting steps
#####################################################################
# Install Creutzml fork
require(devtools)
install_github("creutzml/ffscb", force = T)

# Packages
library(mvtnorm)
library(tidyverse)
library(fda)
library(ffscb)
library(conformalInference.fd)
library(progress)
library(ffscbExtra)
#####################################################################



## Helper functions
#####################################################################
band_score <- function(new_obs, band_lower, band_upper, sig_level) {
  # Calculate the max-width
  max_band_width <- max(band_upper - band_lower, na.rm = T)
  
  # Find any exceedances present
  low_exceedance <- new_obs < band_lower
  high_exceedance <- new_obs > band_upper
  
  # Find the maximum exceedances
  low_max_exceed <- max(((band_lower - new_obs)*low_exceedance)[
    is.finite((band_lower - new_obs)*low_exceedance)
  ])
  high_max_exceed <- max(((new_obs - band_upper)*high_exceedance)[
    is.finite((new_obs - band_upper)*high_exceedance)
  ])
  
  # Calculate the band-score
  band_score <- max_band_width + (2/sig_level)*low_max_exceed + 
    (2/sig_level)*high_max_exceed
  
  # Return the score:
  return(band_score)
}

#####################################################################



## Run a simulation comparing ffscb and conformal
#####################################################################
## Fixed values for the loop
# {err_dist = "t"; n = 50; sigs = 1; dfs = 10;}
# set.seed(101911)

## Parameters to loop over:
n_iters <- 1
err_dist <- "t"
n <- c(30, 100)
dfs <- c(5, 15)
# sigs <- c(1, 3)
sim_parms <- expand.grid(err_dist, n, dfs) #, sigs)
colnames(sim_parms) <- c("err_dist", "n", "dfs") #, "sigs")

# Empty list to store simulation results
sim_results_list <- vector("list", nrow(sim_parms))

# Start the simulation
#####################################################################
# Loop over the simulation parameters
pb <- progress_bar$new(total = nrow(sim_parms)*n_iters)
for (p in 1:nrow(sim_parms)) {
  
  # Empty data frame for results
  sim_results_temp <- data.frame(
    err_dist = rep(sim_parms$err_dist[p], n_iters*4),
    n = rep(sim_parms$n[p], n_iters*4),
    dfs = rep(sim_parms$dfs[p], n_iters*4),
    # sigs = rep(sim_parms$sigs[p], n_iters*4),
    method = vector("character", n_iters*4), 
    prediction = vector("numeric", n_iters*4), 
    out = vector("numeric", n_iters*4),
    out1 = vector("numeric", n_iters*4),
    out2 = vector("numeric", n_iters*4),
    out3 = vector("numeric", n_iters*4),
    # out4 = vector("numeric", n_iters*4),
    # out5 = vector("numeric", n_iters*4),
    band_score = vector("numeric", n_iters*4),
    max_band_width = vector("numeric", n_iters*4)
    # conf_10_out = vector("numeric", n_iters),
    # conf_10_bscore = vector("numeric", n_iters),
    # conf_11_out = vector("numeric", n_iters),
    # conf_11_bscore = vector("numeric", n_iters),
    # ff_10_out = vector("numeric", n_iters),
    # ff_10_bscore = vector("numeric", n_iters),
    # ff_11_out = vector("numeric", n_iters),
    # ff_11_bscore = vector("numeric", n_iters)
  )
  
  # Loop over the iterations
  for (iter in 1:n_iters) {
    
    ## Generate data:
    #################################################################
    # Number of observations and sampling points, 
    n_obs <- sim_parms$n[p]
    n_sp <- 101
    grid <- make_grid(n_sp, rangevals = c(0, 1))
    y_vals_mat <- matrix(0, nrow = n_sp, ncol = n_obs + 2)
    
    # Need an array for the predictor matrix
    x_array <- array(c(matrix(1, nrow = n_sp, ncol = n_obs + 2), 
                       matrix(c(rep(0, n_sp*n_obs/2),
                                rep(1, n_sp*n_obs/2), 
                                rep(0, n_sp), rep(1, n_sp)), 
                              nrow = n_sp, ncol = n_obs + 2)), 
                     dim = c(length(grid), n_obs + 2, 2))
    
    # Beta functions
    b0 <- 1
    b1 <- function(i, t_mat) {
      sin(8*pi*t_mat[i])*exp(-3*t_mat[i]) + t_mat[i]
    }
    
    # True mean functionals
    true_int_mean <- rep(b0, length(grid))
    true_int_beta_mean <- true_int_mean + b1(1:length(grid), grid)
    true_slope_mean <- b1(1:length(grid), grid)
    
    ## Generating response values
    # Random error covariance
    cov.m <- make_cov_m(
      cov.f = ffscb::covf_nonst_matern, #covf_st_matern,
      grid = grid, 
      cov.f.params=c(2, 1/4, 1/4) #c(3/2, 1/4)
    )
    
    # t-distributed errors
    eps_mat <- make_sample(
      mean.v = rep(0, n_sp),
      cov.m = cov.m,
      N = n_obs + 2,
      dist = "rnorm"
    )*sqrt(sim_parms$dfs[p]/rchisq(1, df = sim_parms$dfs[p]))

    
    # Generate the random data
    y_vals_mat <- t(x_array[1,,]%*%t(cbind(true_int_mean, true_slope_mean))) + 
      eps_mat
    
    # matplot(y_vals_mat, type = "l")
    
    
    # Next, remove the last two columns for hold out prediction
    y_vals_mat_temp <- y_vals_mat[,-c(n_obs + 1, n_obs + 2)]
    x_array_temp <- x_array[,-c(n_obs + 1, n_obs + 2),]
    y_vals_mat_ho <- y_vals_mat[,c(n_obs + 1, n_obs + 2)]
    x_array_ho <- x_array[,c(n_obs + 1, n_obs + 2),]
    #################################################################
    
    
    ## Conformal predictions
    #################################################################
    # Before using conformal, we need to convert the forms of data
    grid_list <- list(as.list(grid))
    x_list_temp <- vector("list", length = n_obs)
    y_list_temp <- vector("list", length = n_obs)
    x_list_ho <- vector("list", length = 2)
    
    for (i in 1:n_obs) {
      for (j in 1:n_sp) {
        x_list_temp[[i]][[j]] <- c(x_array_temp[j,i,])
        y_list_temp[[i]][[j]] <- y_vals_mat_temp[j,i]
      }
    }
    
    for (i in 1:2) {
      for (j in 1:n_sp) {
        x_list_ho[[i]][[j]] <- c(x_array_ho[j,i,])
      }
    }
    
    x_list_ho1 <- list(x_list_ho[[1]])
    x_list_ho2 <- list(x_list_ho[[2]])
    
    # Function needed for arguments in conformal
    fun=mean_lists()
    
    # Conformal prediction of c(1,0)
    final.mfData = conformal.fun.split(
      x = x_list_temp, 
      t_x = grid_list,
      y = y_list_temp, 
      t_y = grid_list,
      x0 = x_list_ho1,
      train.fun = fun$train.fun, 
      predict.fun = fun$predict.fun,
      alpha = 0.1,
      split = NULL, 
      seed = FALSE, 
      randomized = FALSE,
      seed.rand = FALSE,
      verbose = FALSE, 
      rho = 0.5,
      s.type = "identity")
    
    # Check coverage
    conf_10_out_temp <- (sum(
      y_vals_mat_ho[,1] < unlist(final.mfData$lo[[1]]) |
        y_vals_mat_ho[,1] > unlist(final.mfData$up[[1]])
    ) > 0)
    
    # Check subinterval coverage
    conf_10_out_temp_by_int <- c()
    for (c in 1:3) {
      start_id <- 34*(c - 1) + 1
      end_id <- 34*c - 1*(c == 3)
      
      conf_10_out_temp_by_int[c] <- (sum(
        y_vals_mat_ho[start_id:end_id,1] < 
          unlist(final.mfData$lo[[1]])[start_id:end_id] |
          y_vals_mat_ho[start_id:end_id,1] > 
          unlist(final.mfData$up[[1]])[start_id:end_id]
      ) > 0)
    }
    
    # Calculate band score
    conf_10_bscore_temp <- band_score(
      new_obs = y_vals_mat_ho[,1], 
      band_lower = unlist(final.mfData$lo[[1]]),
      band_upper = unlist(final.mfData$up[[1]]), 
      sig_level = 0.1
    )
    
    # Calculate mean band width
    conf_10_bwidth_temp <- max(unlist(final.mfData$up[[1]]) - 
                                 unlist(final.mfData$lo[[1]]), 
                               na.rm = T)
    
    # Conformal prediction of c(1,1)
    final.mfData = conformal.fun.split(
      x = x_list_temp, 
      t_x = grid_list,
      y = y_list_temp, 
      t_y = grid_list,
      x0 = x_list_ho2,
      train.fun = fun$train.fun, 
      predict.fun = fun$predict.fun,
      alpha = 0.1,
      split = NULL, 
      seed = FALSE, 
      randomized = FALSE,
      seed.rand = FALSE,
      verbose = FALSE, 
      rho = 0.5,
      s.type = "identity")
    
    # Check coverage
    conf_11_out_temp <- (sum(
      y_vals_mat_ho[,2] < unlist(final.mfData$lo[[1]]) |
        y_vals_mat_ho[,2] > unlist(final.mfData$up[[1]])
    ) > 0)
    
    # Check subinterval coverage
    conf_11_out_temp_by_int <- c()
    for (c in 1:3) {
      start_id <- 34*(c - 1) + 1
      end_id <- 34*c - 1*(c == 3)
      
      conf_11_out_temp_by_int[c] <- (sum(
        y_vals_mat_ho[start_id:end_id,2] < 
          unlist(final.mfData$lo[[1]])[start_id:end_id] |
          y_vals_mat_ho[start_id:end_id,2] > 
          unlist(final.mfData$up[[1]])[start_id:end_id]
      ) > 0)
    }
    
    # Calculate band score
    conf_11_bscore_temp <- band_score(
      new_obs = y_vals_mat_ho[,2], 
      band_lower = unlist(final.mfData$lo[[1]]),
      band_upper = unlist(final.mfData$up[[1]]), 
      sig_level = 0.1
    )
    
    # Calculate mean band width
    conf_11_bwidth_temp <- max(unlist(final.mfData$up[[1]]) - 
                                 unlist(final.mfData$lo[[1]]), 
                               na.rm = T)
    #################################################################
    
    
    ## FFSCB fitting
    #################################################################
    ## And now, if we fit with concurrent regression I made and our bands
    fReg_list <- fRegress_concurrent(y_mat = y_vals_mat_temp, 
                                     x_array = x_array_temp[,,-1])
    
    # First prediction: c(1,0)
    fBands <- predict_concurrent(
      concurrent_list = fReg_list,
      interval = "prediction", 
      err_str = "t",
      new_dat = c(0), 
      conf.level = 0.90, 
      n_int = 3,
      nu0_hat = "singh_df",
      mse_scalar = "ub"
    )
    
    # Check coverage
    ff_10_out_temp <- (sum(
      y_vals_mat_ho[,1] < fBands[[1]][,3] |
        y_vals_mat_ho[,1] > fBands[[1]][,2]
    ) > 0)
    
    # Check subinterval coverage
    ff_10_out_temp_by_int <- c()
    for (c in 1:3) {
      start_id <- 34*(c - 1) + 1
      end_id <- 34*c - 1*(c == 3)
      
      ff_10_out_temp_by_int[c] <- (sum(
        y_vals_mat_ho[start_id:end_id,1] < 
          fBands[[1]][start_id:end_id,3] |
          y_vals_mat_ho[start_id:end_id,1] > 
          fBands[[1]][start_id:end_id,2]
      ) > 0)
    }
    
    # Calculate band score
    ff_10_bscore_temp <- band_score(
      new_obs = y_vals_mat_ho[,1], 
      band_lower = fBands[[1]][,3],
      band_upper = fBands[[1]][,2], 
      sig_level = 0.1
    )
    
    # Calculate mean band width
    ff_10_bwidth_temp <- max(fBands[[1]][,2] - fBands[[1]][,3], 
                             na.rm = T)
    
    # Second prediction: c(1,1)
    fBands2 <- predict_concurrent(
      concurrent_list = fReg_list,
      interval = "prediction", 
      err_str = "t",
      new_dat = c(1), 
      conf.level = 0.90, 
      n_int = 3,
      nu0_hat = "singh_df",
      mse_scalar = "ub"
    )
    
    # Check coverage
    ff_11_out_temp <- (sum(
      y_vals_mat_ho[,2] < fBands2[[1]][,3] |
        y_vals_mat_ho[,2] > fBands2[[1]][,2]
    ) > 0)
    
    # Check subinterval coverage
    ff_11_out_temp_by_int <- c()
    for (c in 1:3) {
      start_id <- 34*(c - 1) + 1
      end_id <- 34*c - 1*(c == 3)
      
      ff_11_out_temp_by_int[c] <- (sum(
        y_vals_mat_ho[start_id:end_id,2] < 
          fBands2[[1]][start_id:end_id,3] |
          y_vals_mat_ho[start_id:end_id,2] > 
          fBands2[[1]][start_id:end_id,2]
      ) > 0)
    }
    
    # Calculate band score
    ff_11_bscore_temp <- band_score(
      new_obs = y_vals_mat_ho[,2], 
      band_lower = fBands2[[1]][,3],
      band_upper = fBands2[[1]][,2], 
      sig_level = 0.1
    )
    
    # Calculate mean band width
    ff_11_bwidth_temp <- max(fBands2[[1]][,2] - fBands2[[1]][,3], 
                             na.rm = T)
    #################################################################
    
    
    ## Update data frame and progress
    #################################################################
    # Update data frame
    sim_results_iter <- matrix(
      c("Conformal", "Conformal", "Fast and Fair", "Fast and Fair",
        "c(1,0)", "c(1,1)", "c(1,0)", "c(1,1)",
        conf_10_out_temp, conf_11_out_temp, ff_10_out_temp, 
        ff_11_out_temp, 
        conf_10_out_temp_by_int[1], conf_11_out_temp_by_int[1], 
        ff_10_out_temp_by_int[1], ff_11_out_temp_by_int[1],
        conf_10_out_temp_by_int[2], conf_11_out_temp_by_int[2], 
        ff_10_out_temp_by_int[2], ff_11_out_temp_by_int[2],
        conf_10_out_temp_by_int[3], conf_11_out_temp_by_int[3], 
        ff_10_out_temp_by_int[3], ff_11_out_temp_by_int[3],
        # conf_10_out_temp_by_int[4], conf_11_out_temp_by_int[4], 
        # ff_10_out_temp_by_int[4], ff_11_out_temp_by_int[4],
        # conf_10_out_temp_by_int[5], conf_11_out_temp_by_int[5], 
        # ff_10_out_temp_by_int[5], ff_11_out_temp_by_int[5],
        conf_10_bscore_temp, conf_11_bscore_temp, ff_10_bscore_temp, 
        ff_11_bscore_temp, 
        conf_10_bwidth_temp, conf_11_bwidth_temp, ff_10_bwidth_temp, 
        ff_11_bwidth_temp), 
      nrow = 4
    )
    
    # Indices for rows
    s_idx <- (iter - 1)*4 + 1
    e_idx <- s_idx + 3
    
    sim_results_temp[s_idx:e_idx, c(4:11)] <- sim_results_iter
    
    # Update progress bar
    pb$tick()
    #################################################################
  }
  
  ## Save data frame to list:
  ###################################################################
  sim_results_list[[p]] <- sim_results_temp
  ###################################################################
}

# Create and save data frame from list
sim_results_df <- do.call(rbind, sim_results_list)
#####################################################################





### This code will reproduce Figure 2 in Creutzinger, Liebl, and 
### Sharp (2024+)
#####################################################################
# Number of observations and sampling points, 
n_obs <- 100
n_sp <- 101
grid <- make_grid(n_sp, rangevals = c(0, 1))
y_vals_mat <- matrix(0, nrow = n_sp, ncol = n_obs + 2)

# Need an array for the predictor matrix
x_array <- array(c(matrix(1, nrow = n_sp, ncol = n_obs + 2), 
                   matrix(c(rep(0, n_sp*n_obs/2),
                            rep(1, n_sp*n_obs/2), 
                            rep(0, n_sp), rep(1, n_sp)), 
                          nrow = n_sp, ncol = n_obs + 2)), 
                 dim = c(length(grid), n_obs + 2, 2))

# Beta functions
b0 <- 1
b1 <- function(i, t_mat) {
  sin(8*pi*t_mat[i])*exp(-3*t_mat[i]) + t_mat[i]
}

# True mean functionals
true_int_mean <- rep(b0, length(grid))
true_int_beta_mean <- true_int_mean + b1(1:length(grid), grid)
true_slope_mean <- b1(1:length(grid), grid)

## Generating response values
# Random error covariance
cov.m <- make_cov_m(
  cov.f = ffscb::covf_nonst_matern, #covf_st_matern,
  grid = grid, 
  cov.f.params=c(2, 1/4, 1/4) #c(3/2, 1/4)
)

# t-distributed errors
eps_mat <- make_sample(
  mean.v = rep(0, n_sp),
  cov.m = cov.m,
  N = n_obs + 2,
  dist = "rnorm"
)*sqrt(15/rchisq(1, df = 15))


# Generate the random data
y_vals_mat_nonst <- t(x_array[1,,]%*%
                        t(cbind(true_int_mean, true_slope_mean))) + 
  eps_mat

# Next, remove the last two columns for hold out prediction
y_vals_mat_temp_nonst <- y_vals_mat_nonst[,-c(n_obs + 1, n_obs + 2)]
x_array_temp_nonst <- x_array[,-c(n_obs + 1, n_obs + 2),]
y_vals_mat_ho_nonst <- y_vals_mat_nonst[,c(n_obs + 1, n_obs + 2)]
x_array_ho_nonst <- x_array[,c(n_obs + 1, n_obs + 2),]

# Random error covariance
cov.m <- make_cov_m(
  cov.f = ffscb::covf_st_matern,
  grid = grid, 
  cov.f.params=c(3/2, 1/4)
)

# t-distributed errors
eps_mat <- make_sample(
  mean.v = rep(0, n_sp),
  cov.m = cov.m,
  N = n_obs + 2,
  dist = "rnorm"
)*sqrt(15/rchisq(1, df = 15))


# Generate the random data
y_vals_mat_st <- t(x_array[1,,]%*%
                        t(cbind(true_int_mean, true_slope_mean))) + 
  eps_mat

# Next, remove the last two columns for hold out prediction
y_vals_mat_temp_st <- y_vals_mat_st[,-c(n_obs + 1, n_obs + 2)]
x_array_temp_st <- x_array[,-c(n_obs + 1, n_obs + 2),]
y_vals_mat_ho_st <- y_vals_mat_st[,c(n_obs + 1, n_obs + 2)]
x_array_ho_st <- x_array[,c(n_obs + 1, n_obs + 2),]


# Quick plot of the curves
y_vals_df_nonst <- y_vals_mat_nonst %>%
  as.data.frame() %>%
  mutate(t = seq(0, 1, length.out = 101), 
         cov_st = "Non-Stationary") %>%
  pivot_longer(-c(t, cov_st), names_to = "Obs", values_to = "Y(t)") %>%
  mutate(`x(t)` = case_when(
    parse_number(Obs) < 51 | parse_number(Obs) == 101 ~ "x(t) = 0",
    parse_number(Obs) > 50 & parse_number(Obs) != 101 ~ "x(t) = 1"
  ))

y_vals_df_st <- y_vals_mat_st %>%
  as.data.frame() %>%
  mutate(t = seq(0, 1, length.out = 101), 
         cov_st = "Stationary") %>%
  pivot_longer(-c(t, cov_st), names_to = "Obs", values_to = "Y(t)") %>%
  mutate(`x(t)` = case_when(
    parse_number(Obs) < 51 | parse_number(Obs) == 101 ~ "x(t) = 0",
    parse_number(Obs) > 50 & parse_number(Obs) != 101 ~ "x(t) = 1"
  ))
y_vals_df <- dplyr::bind_rows(y_vals_df_nonst, y_vals_df_st) %>%
  mutate(cov_st = factor(cov_st, 
                         levels = c("Stationary", "Non-Stationary")))

ggplot() +
  geom_vline(xintercept = c(1/3, 2/3),
             color = "lightgrey") +
  geom_line(aes(x = t, y = `Y(t)`,
                color = `x(t)`, linetype = `x(t)`, group = Obs),
            data = y_vals_df) +
  facet_wrap(vars(cov_st), ncol = 2, 
             strip.position = "top") +
  scale_color_manual(values = c("#D55E00", "#56B4E9")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme_bw(base_size = 16) +
  labs(color = "Predictor:", linetype = "Predictor:") +
  theme(legend.position = c(.5, .9),
        legend.direction = "horizontal",
        legend.background = element_rect(colour = 'black',
                                         fill = 'white',
                                         linetype='solid'),
        # text = element_text(size = 16),
        plot.margin = unit(c(0.1, .7, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text = element_text(size = 16), 
        panel.spacing.x = unit(2.5, "lines")) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = round(c(0, 1/3, 2/3, 1), 2)) +
  coord_cartesian(ylim = c(0, 3))
#####################################################################




### This code will repdroduce Figure 3
### Note: code chunk for Figure 2 needs to be run before this one
#####################################################################
# Quick comparison plot of bands created by conformal vs fast and fair
band_plot_df_y_nonst <- data.frame(
  t = seq(0,1, length.out = 101),
  y_true = y_vals_mat_ho_nonst[,2], 
  cov_st = "Non-Stationary"
)

# Before using conformal, we need to convert the forms of data
grid_list <- list(as.list(grid))
x_list_temp_nonst <- vector("list", length = n_obs)
y_list_temp_nonst <- vector("list", length = n_obs)
x_list_ho_nonst <- vector("list", length = 2)

for (i in 1:n_obs) {
  for (j in 1:n_sp) {
    x_list_temp_nonst[[i]][[j]] <- c(x_array_temp_nonst[j,i,])
    y_list_temp_nonst[[i]][[j]] <- y_vals_mat_temp_nonst[j,i]
  }
}

for (i in 1:2) {
  for (j in 1:n_sp) {
    x_list_ho_nonst[[i]][[j]] <- c(x_array_ho_nonst[j,i,])
  }
}

x_list_ho1_nonst <- list(x_list_ho_nonst[[1]])
x_list_ho2_nonst <- list(x_list_ho_nonst[[2]])

# Function needed for arguments in conformal
fun=mean_lists()

# Conformal prediction of c(1,1)
final.mfData_nonst = conformal.fun.split(
  x = x_list_temp_nonst, 
  t_x = grid_list,
  y = y_list_temp_nonst, 
  t_y = grid_list,
  x0 = x_list_ho2_nonst,
  train.fun = fun$train.fun, 
  predict.fun = fun$predict.fun,
  alpha = 0.1,
  split = NULL, 
  seed = FALSE, 
  randomized = FALSE,
  seed.rand = FALSE,
  verbose = FALSE, 
  rho = 0.5,
  s.type = "identity")

band_plot_df_ci_nonst <- data.frame(
  t = seq(0,1, length.out = 101),
  # y_true = y_vals_mat_ho[,2],
  ci_pred = unlist(final.mfData_nonst$pred),
  ci_lower = unlist(final.mfData_nonst$lo[[1]]),
  ci_upper = unlist(final.mfData_nonst$up[[1]]) #,
  # ff_pred = fBands2[[1]][[1]][,1],
  # ff_lower = fBands2[[1]][[1]][,3],
  # ff_upper = fBands2[[1]][[1]][,2]
) %>%
mutate(band_type = "Conformal \n Inference", 
       cov_st = "Non-Stationary")

## And now, if we fit with concurrent regression I made and our bands
fReg_list_nonst <- fRegress_concurrent(y_mat = y_vals_mat_temp_nonst, 
                                       x_array = x_array_temp_nonst[,,-1])

# prediction: c(1,1)
fBands_nonst <- predict_concurrent(
  concurrent_list = fReg_list_nonst,
  interval = "prediction", 
  err_str = "t",
  new_dat = c(1), 
  conf.level = 0.90, 
  n_int = 3,
  nu0_hat = "singh_df",
  mse_scalar = "ub"
)

band_plot_df_ff_nonst <- data.frame(
  t = seq(0,1, length.out = 101),
  # y_true = y_vals_mat_ho[,2],
  # pred = unlist(final.mfData$pred),
  # lower = unlist(final.mfData$lo[[1]]),
  # upper = unlist(final.mfData$up[[1]]) #,
  ff_pred = fBands_nonst[[1]][,1],
  ff_lower = fBands_nonst[[1]][,3],
  ff_upper = fBands_nonst[[1]][,2]
) %>%
  mutate(band_type = "Fast and \n Fair", 
         cov_st = "Non-Stationary")

# Quick comparison plot of bands created by conformal vs fast and fair
band_plot_df_y_st <- data.frame(
  t = seq(0,1, length.out = 101),
  y_true = y_vals_mat_ho_st[,2], 
  cov_st = "Stationary"
)

# Before using conformal, we need to convert the forms of data
grid_list <- list(as.list(grid))
x_list_temp_st <- vector("list", length = n_obs)
y_list_temp_st <- vector("list", length = n_obs)
x_list_ho_st <- vector("list", length = 2)

for (i in 1:n_obs) {
  for (j in 1:n_sp) {
    x_list_temp_st[[i]][[j]] <- c(x_array_temp_st[j,i,])
    y_list_temp_st[[i]][[j]] <- y_vals_mat_temp_st[j,i]
  }
}

for (i in 1:2) {
  for (j in 1:n_sp) {
    x_list_ho_st[[i]][[j]] <- c(x_array_ho_st[j,i,])
  }
}

x_list_ho1_st <- list(x_list_ho_st[[1]])
x_list_ho2_st <- list(x_list_ho_st[[2]])

# Function needed for arguments in conformal
fun=mean_lists()

# Conformal prediction of c(1,1)
final.mfData_st = conformal.fun.split(
  x = x_list_temp_st, 
  t_x = grid_list,
  y = y_list_temp_st, 
  t_y = grid_list,
  x0 = x_list_ho2_st,
  train.fun = fun$train.fun, 
  predict.fun = fun$predict.fun,
  alpha = 0.1,
  split = NULL, 
  seed = FALSE, 
  randomized = FALSE,
  seed.rand = FALSE,
  verbose = FALSE, 
  rho = 0.5,
  s.type = "identity")

band_plot_df_ci_st <- data.frame(
  t = seq(0,1, length.out = 101),
  # y_true = y_vals_mat_ho[,2],
  ci_pred = unlist(final.mfData_st$pred),
  ci_lower = unlist(final.mfData_st$lo[[1]]),
  ci_upper = unlist(final.mfData_st$up[[1]]) #,
  # ff_pred = fBands2[[1]][[1]][,1],
  # ff_lower = fBands2[[1]][[1]][,3],
  # ff_upper = fBands2[[1]][[1]][,2]
) %>%
  mutate(band_type = "Conformal \n Inference", 
         cov_st = "Stationary")

## And now, if we fit with concurrent regression I made and our bands
fReg_list_st <- fRegress_concurrent(y_mat = y_vals_mat_temp_st, 
                                       x_array = x_array_temp_st[,,-1])

# prediction: c(1,1)
fBands_st <- predict_concurrent(
  concurrent_list = fReg_list_st,
  interval = "prediction", 
  err_str = "t",
  new_dat = c(1), 
  conf.level = 0.90, 
  n_int = 3,
  nu0_hat = "singh_df",
  mse_scalar = "ub"
)

band_plot_df_ff_st <- data.frame(
  t = seq(0,1, length.out = 101),
  # y_true = y_vals_mat_ho[,2],
  # pred = unlist(final.mfData$pred),
  # lower = unlist(final.mfData$lo[[1]]),
  # upper = unlist(final.mfData$up[[1]]) #,
  ff_pred = fBands_st[[1]][,1],
  ff_lower = fBands_st[[1]][,3],
  ff_upper = fBands_st[[1]][,2]
) %>%
  mutate(band_type = "Fast and \n Fair", 
         cov_st = "Stationary")

## Combine data frames
band_plot_df_ci <- dplyr::bind_rows(band_plot_df_ci_nonst,
                                    band_plot_df_ci_st) %>%
  mutate(cov_st = factor(cov_st, 
                         levels = c("Stationary", "Non-Stationary")))
band_plot_df_ff <- dplyr::bind_rows(band_plot_df_ff_nonst, 
                                    band_plot_df_ff_st) %>%
  mutate(cov_st = factor(cov_st, 
                         levels = c("Stationary", "Non-Stationary")))
band_plot_df_y <- dplyr::bind_rows(band_plot_df_y_nonst, 
                                   band_plot_df_y_st) %>%
  mutate(cov_st = factor(cov_st, 
                         levels = c("Stationary", "Non-Stationary")))

# # Plot the comparison
fct_color <- c("black", "#D55E00", "#56B4E9")
fct_linetype <- c(1, 2, 4)
ggplot() +
  geom_vline(xintercept = c(1/3, 2/3),
             color = "lightgrey") +
  geom_ribbon(aes(x = t, ymin = ff_lower, ymax = ff_upper,
                  fill = band_type, linetype = band_type),
              alpha = 0.5, linewidth = 1.25, #linetype = "solid",
              data = band_plot_df_ff, key_glyph = "blank") +
  geom_ribbon(aes(x = t, ymin = ci_lower, ymax = ci_upper,
                  color = band_type, linetype = band_type),
              # linetype = "dotdash",
              fill = NA,
              alpha = 0.5,
              # linetype = "dotdash",
              linewidth = 1.25,
              data = band_plot_df_ci, key_glyph = "blank") +
  geom_line(aes(x = t, y = y_true, group = cov_st),
            color = "black", linewidth = 1.25,
            data = band_plot_df_y) +
  geom_line(aes(x = t, y = ff_pred, color = band_type, 
                linetype = band_type),
            linewidth = 1.25,
            data = band_plot_df_ff) +
  facet_wrap(vars(cov_st), nrow = 2, 
             strip.position = "right") +
  theme_bw(base_size = 20) +
  # theme(text = element_text(size = 16)) +
  labs(#color = "Prediction Band: ",
    # linetype = "Prediction Band: ",
    # fill = "Prediction Band: ",
    y = "Y(t)") +
  scale_color_manual(name = "Prediction Band: ",
                     breaks = c("Conformal \n Inference",
                                "Fast and \n Fair"),
                     values = c("Conformal \n Inference" = "#D55E00",
                                "Fast and \n Fair" = "#56B4E9")) +
  scale_fill_manual(name = "Prediction Band: ",
                    breaks = c("Conformal \n Inference",
                               "Fast and \n Fair"),
                    values = c("Conformal \n Inference" = NA,
                               "Fast and \n Fair" = "#56B4E9")) +
  scale_linetype_manual(name = "Prediction Band: ",
                        breaks = c("Conformal \n Inference",
                                   "Fast and \n Fair"),
                        values = c("Conformal \n Inference" = "dotdash",
                                   "Fast and \n Fair" = "solid")) +
  guides(fill = "none",
         color = guide_legend(override.aes = list(fill = NA))) +
  # title = "Comparison of 90% Prediction Bands for x(t) = 1",
  # subtitle = expression(nu[0]*" = 30, n = 100")) +
  theme(legend.position = c(.5, .95),
        legend.direction = "horizontal",
        legend.background = element_rect(colour = 'black',
                                         fill = 'white',
                                         linetype='solid'),
        # text = element_text(size = 16),
        plot.margin = unit(c(1, 1, 0.1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.y = element_text(size = 20)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = round(c(0, 1/3, 2/3, 1), 2)) +
  coord_cartesian(ylim = c(-0.25, 3))
#####################################################################

