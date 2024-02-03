#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### May 23rd, 2023                                                ###
###                                                               ###
###   Creating simultaneous prediction bands for a new            ###
### observation. Specific interest will be for those in the       ###
### amputee group. Model will be built with best predictors as    ###
### well, and those will be predicted at thier average levels,    ###
### aside from categorical variables. Those will be estimated     ###
### separately when possible.                                     ###
#####################################################################


## Creutzinger version of package
require(devtools)
install_github("creutzml/ffscbExtra")

## Packages
library(tidyverse)
library(ffscbExtra)
library(glmnet)
require(methods)
library(abind)
library(conformalInference.fd)
library(gridExtra)
library(cowplot)




## Initial setup variables
#####################################################################
grid <- seq(from = 0, to = 1, len=101)
time_vec <- grid

## Load data in:
dir_path <- file.path(here::here())
## 1/28/2024: need to wait to see if it's okay to upload data to repo

## Choose one of different functional outcome variables:
slct <- 7 # <- focusing on "FRONT_V" looks most promising to me
Y_curves <- readr::read_delim(
  file = c("data/BACK_AP.csv", # 1
           "data/BACK_ML.csv", # 2
           "data/BACK_V.csv",  # 3
           "data/BACK_RES.csv",# 4
           "data/FRONT_AP.csv",# 5
           "data/FRONT_ML.csv",# 6
           "data/FRONT_V.csv", # 7
           "data/FRONT_RES.csv"# 8
  )[slct], 
  delim=";",
  col_names = FALSE) %>% 
  as.matrix(.)


## Classic predictors
X_pred <- readr::read_delim(
  file = "DATA/Covariates.csv", 
  delim=";",
  col_names = TRUE) 
X_pred <- X_pred[,-1]
colnames(X_pred) <- gsub("'", "", colnames(X_pred))
X_pred$Sex <- factor(X_pred$Sex, 
                     levels = c(-1, 1),
                     labels = c("female", "male"))
# All amputees are male; make sure this is reflected in data
X_pred$Sex[155:161] <- "male"

### May 24th, 2023: there was confusion on the amputee sprinters, 
### because the data shows that two of them are female, but in 
### Willwacher et al. (2016), they are described to all be male. When
### Dominik asked Steffen, Steffen was confident this was just a typo
### in the data. So I will change the Sex of the two "female" amp's 
### to male. 
# # Need to find the indices of amputee sprinters
# amp_rel_pr <- c(105.9, 101.2, 102.4, 116.9, 112.4, 110.3, 116.1)
# amp_indices <- with(X_pred, {
#   which(`Relative PR` %in% amp_rel_pr)
# })
# View(X_pred[amp_indices,])


# Add a predictor for amputee or not
N <- n_obs <- ncol(Y_curves)
X_pred$Amp <- factor(c(rep(0, N - 7), rep(1, 7)),
                     levels = c(0, 1),
                     labels = c("non-amputee", "amputee"))
#####################################################################



### Need to run phase alignment, based on the second peak of the
### force curve
#####################################################################
## Start by finding the optimal point of all curves between 0.7, 0.9
int_start <- round(0.7*101)
int_end <- round(0.9*101)
max_idx <- apply(Y_curves[(int_start:int_end),], 
                 MARGIN = 2,
                 FUN = which.max) + int_start - 1
max_vals <- vector("numeric", n_obs)
for (i in 1:n_obs) {
  max_vals[i] <- Y_curves[max_idx[i],i]
}
max_pts <- cbind((max_idx - 1)/101, max_vals) %>% as.data.frame()
# matplot(x = 1:101,  y = Y_curves, type="l", lty=1, col="gray")
# points(max_pts[,1], max_pts[,2], col = "red")

# Long form data for plotting
Y_curves_dat <- Y_curves %>%
  as.data.frame() %>%
  dplyr::mutate(t_pts = grid) %>%
  tidyr::pivot_longer(cols = -t_pts, 
                      names_to = "Observation", 
                      values_to = "Front_V")

## Recreates Figure 7 in Chapter 3, Creutzinger (2024) dissertation
# Plot of original Front V curves, with a point highlighting their max
y_max_plot <- ggplot() + 
  geom_line(aes(x = 100*t_pts, y = Front_V, group = Observation), 
            color = "gray", 
            data = Y_curves_dat) +
  geom_point(aes(x = 100*V1, y = max_vals), 
             data = max_pts, color = "red") +
  scale_x_continuous(breaks = round(seq(0, 100, 100/3), 2), 
                     expand = c(0,0)) +
  labs(x = "% Push-Off Phase", 
       y = "Front Vertical Force (N/kg)") +
  theme_bw() +
  theme(text = element_text(size = 16),
        plot.margin = unit(c(0.6, 1, .6, .6), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


## How does a plot look if we align at those points? 
Y_curves_shifted <- Y_curves %>%
  as.data.frame() %>%
  dplyr::mutate(t_pts = 1:101) %>%
  tidyr::pivot_longer(cols = -t_pts, 
                      names_to = "Obs", 
                      values_to = "Front_V")

## Phase shift
for (i in 1:n_obs) {
  Y_curves_shifted[Y_curves_shifted$Obs == paste0("X", i),]$t_pts <-
    Y_curves_shifted[Y_curves_shifted$Obs == paste0("X", i),]$t_pts - 
    max_idx[i]
}

# Plot of how it looks
ggplot() +
  geom_line(aes(x = t_pts, y = Front_V, group = Obs),
            data = Y_curves_shifted)

## Convert back to wide format, then interpolate starting values
# Need to have the order of columns
col_names_order <- paste0("X", 1:161)
# Wide format
Y_curves_shifted_wide <- Y_curves_shifted %>%
  dplyr::filter(t_pts <= 0) %>%
  tidyr::pivot_wider(names_from = "Obs", values_from = "Front_V") %>%
  dplyr::arrange(t_pts)
  # dplyr::mutate(t_pts = 1:74) %>%
  # dplyr::select(t_pts, all_of(col_names_order))

# Interpolate values to get 101 sampling points for all observations
Y_curves_interp <- apply(
  Y_curves_shifted_wide[,-1], 2, FUN = function(x) {
    approx(x, n = 101)$y
  })

# Long format
Y_curves_interp_long <- Y_curves_interp %>%
  as.data.frame() %>%
  dplyr::mutate(t_pts = grid) %>%
  tidyr::pivot_longer(cols = -t_pts, 
                      names_to = "Obs", 
                      values_to = "Front_V")

# matplot(Y_curves_shifted_wide[,-1], type = "l")
y_shifted_plot <- ggplot() + 
  geom_line(aes(x = 100*t_pts, y = Front_V, group = Obs), 
            data = Y_curves_interp_long, 
            color = "gray") + 
  scale_x_continuous(breaks = round(seq(0, 100, 100/3), 2), 
                     expand = c(0,0)) +
  labs(x = "% Push-Off Phase", 
       y = "Force (N/kg)") +
  theme_bw() +
  theme(text = element_text(size = 16), 
        plot.margin = unit(c(0.6, 1, .6, .6), "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# ggsave(filename = "original_responses.pdf",
#        plot = y_max_plot,
#        device = "pdf",
#        path = paste0("/Users/creutzml/Library/Mobile Documents/com",
#                      "~apple~CloudDocs/Documents/Dissertation/func",
#                      "tional_data_analysis/data/Sprint_Start/figures"),
#        width = 11,
#        height = 8.5,
#        units = "in")

plot_grid(y_max_plot, y_shifted_plot, ncol = 1, align = "v")
#####################################################################



### Plot a tau estimate for the phase aligned Front_V force
#####################################################################
# Create matrix from shifted curves and new grid
Y_curves_shifted_mat <- as.matrix(Y_curves_interp)
n_sp <- nrow(Y_curves_shifted_mat)
grid <- seq(0, 1, length.out = n_sp)

# Estimate and plot tau
front_v_tau_est <- ffscb::tau_fun(Y_curves_shifted_mat)
ggplot() +
  geom_line(aes(x = grid, y = front_v_tau_est)) +
  theme_bw() + 
  labs(title = expression(tau*" Estimate for Front V Kinetic"), 
       x = "T", 
       y = expression(tau)) +
  theme(text = element_text(size=16))
#####################################################################




### Exploring the predictor variables
#####################################################################
## We choose a point along the grid to isolate a random set
## of y values, upon which we can run classic lasso over all the 
## covariates
# Create the data structure needed
sp_iso <- seq(1, n_sp, by = 2)
Y_vals <- Y_curves_shifted_mat[sp_iso,]
X_pred_mat <- as.matrix(X_pred)

var_selection <- as.data.frame(matrix(
  0, nrow = length(sp_iso), ncol = ncol(X_pred) + 1
))
colnames(var_selection) <- c("sp_iso", colnames(X_pred))
var_selection$sp_iso <- sp_iso

for (i in 1:length(sp_iso)) {
  # Run the lasso
  cv_model <- cv.glmnet(
    x = data.matrix(X_pred), y = c(Y_vals[i,]), alpha = 1
  )
  se1_lambda <- cv_model$lambda.1se

  # Grab the results
  se1_model <- glmnet(data.matrix(X_pred), c(Y_vals[i,]),
                      alpha = 1, lambda = se1_lambda)
  test_se1 <- coef(se1_model)

  # Mark down the variables chosen
  var_names <- names(test_se1[test_se1[,1] != 0,])
  col_idx <- which(colnames(var_selection) %in% var_names)
  var_selection[i, col_idx] <- var_selection[i, col_idx] + 1
}

# See what the results are
sort(colSums(var_selection), decreasing = T)
#####################################################################




### Fit model with amputee, demographics, and TPush_V
#####################################################################
# Create necessary model objects
n_obs <- ncol(Y_curves_shifted_mat)
n_sp <- nrow(Y_curves_shifted_mat)
resp_mat <- Y_curves_shifted_mat[,which(X_pred$Amp == "non-amputee")]
n_obs_nonamp <- ncol(resp_mat)
X_pred_nonamp <- X_pred[which(X_pred$Amp == "non-amputee"),]

## Prep the factor variables
# Amputee: amputees = 1, non-amputees = 0
amp_var <- as.numeric(X_pred$Amp)
amp_var[amp_var == 1] <- 0
amp_var[amp_var == 2] <- 1

# Sex: males = 1, females = 0
sex_var <- as.numeric(X_pred_nonamp$Sex)
sex_var[sex_var == 1] <- 0
sex_var[sex_var == 2] <- 1

pred_array <- abind(
  # matrix(rep(as.numeric(amp_var), n_sp), 
  #        ncol = n_obs, 
  #        nrow = n_sp, 
  #        byrow = T), 
  matrix(rep(as.numeric(X_pred_nonamp$Mass), n_sp), 
         ncol = n_obs_nonamp, 
         nrow = n_sp, 
         byrow = T), 
  matrix(rep(as.numeric(X_pred_nonamp$Age), n_sp),
         ncol = n_obs_nonamp,
         nrow = n_sp,
         byrow = T),
  matrix(rep(as.numeric(X_pred_nonamp$Height), n_sp),
         ncol = n_obs_nonamp,
         nrow = n_sp,
         byrow = T),
  matrix(rep(as.numeric(sex_var), n_sp),
         ncol = n_obs_nonamp,
         nrow = n_sp,
         byrow = T),
  matrix(rep(as.numeric(X_pred_nonamp$TPush_V), n_sp), 
         ncol = n_obs_nonamp, 
         nrow = n_sp, 
         byrow = T),
  along = 3
)

# Fit the concurrent regression model
nonamp_model <- fRegress_concurrent(resp_mat, pred_array)

front_v_tau_est <- ffscb::tau_fun(nonamp_model$y_resid_mat)
tau_realigned_plot <- ggplot() +
  geom_line(aes(x = 100*grid, y = front_v_tau_est)) +
  scale_x_continuous(breaks = round(seq(0, 100, 100/3), 2), 
                     expand = c(0,0)) +
  theme_bw() + 
  labs(y = expression("Roughness estimate "*hat(tau)), 
       x = "% Push-Off Phase") +
  theme(text = element_text(size=16), 
        plot.margin = unit(c(1, 1, 1, .5), "cm"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


## September 12th, 2023: fit the same model for Conformal Predictions
# Need to create the empty lists necessary to use conformal
resp_list <- vector("list", length = ncol(resp_mat))
design_list <- vector("list", length = ncol(resp_mat))

for (i in 1:ncol(resp_mat)) {
  for (j in 1:nrow(resp_mat)) {
    design_list[[i]][[j]] <- c(pred_array[j,i,])
    resp_list[[i]][[j]] <- resp_mat[j,i]
  }
}

grid_list <- list(as.list(grid))

# Function needed for arguments in conformal
fun=mean_lists()

# Holdout design matrix values
X_pred_amp <- X_pred[which(X_pred$Amp == "amputee"), ]
X_pred_amp_model <- X_pred_amp[,c(2,3,4,6,10)]
design_list_ho <- vector("list", length = nrow(X_pred_amp_model))

for (i in 1:nrow(X_pred_amp_model)) {
  for (j in 1:nrow(resp_mat)) {
    design_list_ho[[i]][[j]] <- as.numeric(X_pred_amp_model[i,])
  }
}

# Fit the conformal model and get predictions
nonamp_conformal_model1 <- conformal.fun.split(
  x = design_list, 
  t_x = grid_list,
  y = resp_list, 
  t_y = grid_list,
  x0 = list(design_list_ho[[1]]),
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

nonamp_conformal_model2 <- conformal.fun.split(
  x = design_list, 
  t_x = grid_list,
  y = resp_list, 
  t_y = grid_list,
  x0 = list(design_list_ho[[2]]),
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

nonamp_conformal_model3 <- conformal.fun.split(
  x = design_list, 
  t_x = grid_list,
  y = resp_list, 
  t_y = grid_list,
  x0 = list(design_list_ho[[3]]),
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

nonamp_conformal_model4 <- conformal.fun.split(
  x = design_list, 
  t_x = grid_list,
  y = resp_list, 
  t_y = grid_list,
  x0 = list(design_list_ho[[4]]),
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

nonamp_conformal_model5 <- conformal.fun.split(
  x = design_list, 
  t_x = grid_list,
  y = resp_list, 
  t_y = grid_list,
  x0 = list(design_list_ho[[5]]),
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

nonamp_conformal_model6 <- conformal.fun.split(
  x = design_list, 
  t_x = grid_list,
  y = resp_list, 
  t_y = grid_list,
  x0 = list(design_list_ho[[6]]),
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

nonamp_conformal_model7 <- conformal.fun.split(
  x = design_list, 
  t_x = grid_list,
  y = resp_list, 
  t_y = grid_list,
  x0 = list(design_list_ho[[7]]),
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

# Empty vectors to store results
ff_out_vec <- vector("logical", length = 7)
ci_out_vec <- vector("logical", length = 7)

#####################################################################



### Quick look at the functional residuals and estimated coefficients
#####################################################################
## Recreates Figure 9 in Chapter 3, Creutzinger (2024) dissertation
# MVN test of the residuals:
library(MVN)
MVN::mvn(nonamp_model$y_resid_mat, mvnTest = "mardia", 
         multivariatePlot = "qq")

# Plot of residuals for curiosity sake:
matplot(nonamp_model$y_resid_mat, nonamp_model$y_hat_mat, 
        type = "l",
        main = "Plot of Functional Residuals vs Fitted Values", 
        xlab = "Residuals", 
        ylab = "Fitted Values")

## How do the confidence bands for the beta parameters look?
beta_bands <- confint_concurrent(nonamp_model, 
                                 n_int = 3, 
                                 conf.level = .90)
## Notes: 
##  - Age is not significant anywhere along the domain
##  - Height is not significant anywhere along the domain
##  - Sex is not significant anywhere along the domain

## Recreates Figure 8 in Chapter 3, Creutzinger (2024) dissertation
# Plot the estimated coefficient bands
par(mfrow = c(2, 3), 
    mar = c(3,3,3,3))
# plot(beta_bands[[2]],
#      main = paste0("95% FFSCB for Estimated \n",
#                    " Coefficient of Amputee"))
# abline(h = 0, lty = "dashed")
plot(beta_bands[[2]], 
     main = paste0("90% FFSCB for Estimated \n",
                   " Coefficient of Mass"))
abline(h = 0, lty = "dashed")
plot(beta_bands[[3]], 
     main = paste0("90% FFSCB for Estimated \n",
                   " Coefficient of Age"))
abline(h = 0, lty = "dashed")
plot(beta_bands[[4]], 
     main = paste0("90% FFSCB for Estimated \n",
                   " Coefficient of Height"))
abline(h = 0, lty = "dashed")
plot(beta_bands[[5]], 
     main = paste0("90% FFSCB for Estimated \n",
                   " Coefficient of Sex"))
abline(h = 0, lty = "dashed")
plot(beta_bands[[6]], 
     main = paste0("90% FFSCB for Estimated \n",
                   " Coefficient of TPush_V"))
abline(h = 0, lty = "dashed")
#####################################################################



### Now, we try to predict a sprinter's Front\_V force with the 
### demographics of each amputee sprinter
#####################################################################
# What if we predict the non-amputee counterpart of each individual
# amputee?
par(mfrow = c(2, 4))
X_pred_amp <- X_pred[which(X_pred$Amp == "amputee"), ]
X_pred_amp[, c(2,3,4,6,10)]

# Before we predict the bands, let's use a distance matrix calculated
# with Mahalanobis distance to find the nearest neighbor of each 
# amputee sprinter
cholMaha <- function(X) {
  dec <- chol( cov(X) )
  tmp <- forwardsolve(t(dec), t(X) )
  dist(t(tmp))
}
X_pred_md <- cholMaha(X_pred[, c(2,3,4,10)])
X_pred_amp_md <- as.matrix(X_pred_md)[-c(155:161),155:161]

# Find the index for the observation with minimum distance to each
X_pred_nonamp_sim <- apply(X_pred_amp_md, 2, FUN = which.min)

## First amputee runner
X_pred_amp1 <- X_pred_amp[1, c(2,3,4,6,10)]

# First amputee
new_dat_mat1 <- matrix(
  rep(c(X_pred_amp1$Mass[1], 
        X_pred_amp1$Age[1], 
        X_pred_amp1$Height[1], 
        1, 
        X_pred_amp1$TPush_V[1]), 
      each = n_sp), 
  nrow = n_sp
)
non_amp_pred_band1 <- predict_concurrent(
  concurrent_list = nonamp_model, 
  interval = "prediction", 
  new_dat = new_dat_mat1, 
  err_str = "t", 
  conf.level = 0.90,
  n_int = 3, 
  nu0_hat = "singh_df",
  mse_scalar = "ub"
)

### Plot of the mean and band
# Get the band
non_amp_pred_band_data1 <- non_amp_pred_band1[[1]]
non_amp_pred_band_data1_plot <- data.frame(
  t_pts = grid,
  mean = non_amp_pred_band_data1[,1],
  upper = non_amp_pred_band_data1[,2],
  lower = non_amp_pred_band_data1[,3],
  amp = Y_curves_shifted_mat[,155], 
  nonamp = resp_mat[,X_pred_nonamp_sim[1]]
) %>%
  pivot_longer(cols = -t_pts, 
               values_to = "Front_V", 
               names_to = "Estimate")

# Set the colors
est_colors <- c("upper" = "#E69F00", 
                "lower" = "#E69F00",
                "mean" = "#E69F00", 
                "amp" = "#56B4E9", 
                "nonamp" = "black")
est_lines <- c("upper" = "dashed", 
               "lower" = "dashed",
               "mean" = "solid", 
               "amp" = "solid", 
               "nonamp" = "solid")

amp_plot1 <- ggplot() +
  geom_line(aes(x = 100*t_pts, y = Front_V, 
                color = Estimate, linetype = Estimate), 
            data = non_amp_pred_band_data1_plot) +
  scale_color_manual(values = est_colors) +
  scale_linetype_manual(values = est_lines) +
  scale_x_continuous(breaks = seq(0,100,20)) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 20)) +
  labs(x = "% Push-Off Phase (realigned)", 
       y = "Front Vertical Force (N/kg)", 
       title = "First Amputee") +
  coord_cartesian(ylim = c(0, 1250))
# title = paste0("95% Simultaneous Prediction Band to Predict ",
#                "1st Amputee"),
# subtitle = paste0("Estimated for a ",
#                   "sprinter with Mass = 85.7kg, Age = 35 ", 
#                   "years, Height = 2m, Sex = Male, and ",
#                   "TPush_V = 0.441"))

# Plot of "fair" critical value
# non_amp_pred_band_crit1 <- non_amp_pred_band1[[2]] %>%
#   data.frame() %>%
#   dplyr::rename("u_hat" = ".") %>%
#   dplyr::mutate(t_pts = grid)
# 
# uhat_plot1 <- ggplot() +
#   geom_vline(xintercept = c(33 + 1/3, 66 + 2/3), 
#              color = "lightgray") +
#   geom_line(aes(x = 100*t_pts, y = u_hat), 
#             data = non_amp_pred_band_crit1) +
#   scale_x_continuous(breaks = round(seq(0, 100, 33 + 1/3), 2), 
#                      expand = c(0,0)) +
#   theme_bw() + 
#   theme(text = element_text(size = 16), 
#         plot.margin = unit(c(.5, 1, 1, 1), "cm"), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) +
#   labs(x = "% Push-Off Phase",
#        y = expression("Critical value "*hat(u)[t[nu[0]]*", "*alpha/2]^"*"))
# 
# library(cowplot)
# plot_grid(tau_realigned_plot, uhat_plot1, 
#           align = "v", nrow = 2, labels = c("(a)", "(b)"))

### September 12th, 2023
## Creating a comparison plot of Conformal Inference vs Fast and Fair
non_amp_pred_band_data1_comp_plot <- data.frame(
  t_pts = grid,
  ff_mean = non_amp_pred_band_data1[,1],
  ff_upper = non_amp_pred_band_data1[,2],
  ff_lower = non_amp_pred_band_data1[,3],
  ci_mean = unlist(nonamp_conformal_model1$pred[[1]]),
  ci_upper = unlist(nonamp_conformal_model1$up[[1]]),
  ci_lower = unlist(nonamp_conformal_model1$lo[[1]]),
  amp = Y_curves_shifted_mat[,155], 
  nonamp = resp_mat[,X_pred_nonamp_sim[1]]
) 

# Does the amputee exceed the bands?
ff_out_vec[1] <- with(non_amp_pred_band_data1_comp_plot, {
  sum(amp < ff_lower | amp > ff_upper)
})
ci_out_vec[1] <- with(non_amp_pred_band_data1_comp_plot, {
  sum(amp < ci_lower | amp > ci_upper)
})

# Comparison plot
conf_ff_comp1 <- ggplot() +
  geom_vline(xintercept = c(33 + 1/3, 66 + 2/3), color = "lightgray") +
  geom_ribbon(aes(x = 100*t_pts, ymin = ff_lower, ymax = ff_upper), 
              fill = "#56B4E9", alpha = 0.5, linewidth = 1.25,
              data = non_amp_pred_band_data1_comp_plot) +
  geom_ribbon(aes(x = 100*t_pts, ymin = ci_lower, ymax = ci_upper), 
              fill = NA, alpha = 0.5, color = "#D55E00",
              linewidth = 1.25, linetype = "dotdash",
              data = non_amp_pred_band_data1_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = nonamp), 
            color = "black", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data1_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = amp), 
            color = "black", linetype = "dashed", linewidth = 1.25,
            data = non_amp_pred_band_data1_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = ff_mean), 
            color = "#56B4E9", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data1_comp_plot) +
  scale_x_continuous(breaks = round(seq(0,100,(33 + 1/3)), 2)) +
  coord_cartesian(ylim = c(0, 1200), expand = FALSE) +
  theme_bw() +
  labs(x = "% Push-Off Phase", 
       y = "Force (N/kg)") + #, 
  # title = "First Amputee") +
  theme( #legend.position = "none", 
    text = element_text(size = 20), 
    plot.margin = unit(c(.5, .9, .5, .5), "cm"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  annotate("label", 
           x = 35, y = 1100, 
           label = "First Amputee", 
           size = 6) 

# # Mean
# par(mar = c(3,3,3,3))
# plot(grid, non_amp_pred_band_data1[,1],
#      type = "l", xlab = "T", ylab = "Front V Force",
#      xlim = c(0, .90),
#      ylim = c(0, max(Y_curves + 20)), 
#      main = paste0("95% Simultaneous Prediction Band for 1st Amputee"))
# # Non-amputees
# for (i in nonamp_idx1) {
#   lines(x = grid,
#         y = resp_mat[, i],
#         col = "gray")
# }
# # Band
# lines(x = grid, y = non_amp_pred_band_data1[,1],
#       col = "black", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data1[,2],
#       col = "red", lty = "dashed", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data1[,3],
#       col = "red", lty = "dashed", lwd = 2)
# lines(x = grid, y = Y_curves_shifted_mat[,155], col = "green", lwd = 3)


## Second amputee runner
X_pred_amp2 <- X_pred_amp[2, c(2,3,4,6,10)]

# Non-amputee runners that are at least 30 years old and male
nonamp_idx2 <- which(X_pred_nonamp$Mass < 75 &
                       X_pred_nonamp$Mass > 70 &
                       X_pred_nonamp$Sex == "male" &
                       X_pred_nonamp$TPush_V > 0.4)
# First amputee
new_dat_mat2 <- matrix(
  rep(c(X_pred_amp2$Mass[1], 
        X_pred_amp2$Age[1], 
        X_pred_amp2$Height[1], 
        1, 
        X_pred_amp2$TPush_V[1]), 
      each = n_sp), 
  nrow = n_sp
)
non_amp_pred_band2 <- predict_concurrent(
  concurrent_list = nonamp_model, 
  interval = "prediction", 
  new_dat = new_dat_mat2, 
  err_str = "t", 
  conf.level = 0.90,
  n_int = 3, 
  nu0_hat = "singh_df",
  mse_scalar = "ub"
)

### Plot of the mean and band
# Get the band
non_amp_pred_band_data2 <- non_amp_pred_band2[[1]][[1]]
non_amp_pred_band_data2_plot <- data.frame(
  t_pts = grid,
  mean = non_amp_pred_band_data2[,1],
  upper = non_amp_pred_band_data2[,2],
  lower = non_amp_pred_band_data2[,3],
  amp = Y_curves_shifted_mat[,156], 
  nonamp = resp_mat[,X_pred_nonamp_sim[2]]
) %>%
  pivot_longer(cols = -t_pts, 
               values_to = "Front_V", 
               names_to = "Estimate") 

amp_plot2 <- ggplot() +
  geom_line(aes(x = 100*t_pts, y = Front_V, 
                color = Estimate, linetype = Estimate), 
            data = non_amp_pred_band_data2_plot) +
  scale_color_manual(values = est_colors) +
  scale_linetype_manual(values = est_lines) +
  scale_x_continuous(breaks = seq(0,100,20)) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 20)) +
  labs(x = "% Push-Off Phase (realigned)", 
       y = "Front Vertical Force (N/kg)",
       title = "Second Amputee") +
  coord_cartesian(ylim = c(0, 1250))
# title = paste0("95% Simultaneous Prediction Band to Predict ",
#                "2nd Amputee"),
# subtitle = paste0("Estimated for a ",
#                   "sprinter with Mass = 73.8 kg, Age = 32 ", 
#                   "years, Height = 1.89 m, Sex = Male, and ",
#                   "TPush_V = 0.475"))

# Plot of "fair" critical value
non_amp_pred_band_crit2 <- non_amp_pred_band2[[1]][[2]] %>%
  data.frame() %>%
  dplyr::rename("u_hat" = ".") %>%
  dplyr::mutate(t_pts = grid)

uhat_plot2 <- ggplot() +
  geom_line(aes(x = t_pts, y = u_hat), 
            data = non_amp_pred_band_crit2) +
  theme_bw() + 
  labs(x = "T", y = expression(hat(u)[alpha/2]), 
       title = "Second Amputee")

### September 12th, 2023
## Creating a comparison plot of Conformal Inference vs Fast and Fair
non_amp_pred_band_data2_comp_plot <- data.frame(
  t_pts = grid,
  ff_mean = non_amp_pred_band_data2[,1],
  ff_upper = non_amp_pred_band_data2[,2],
  ff_lower = non_amp_pred_band_data2[,3],
  ci_mean = unlist(nonamp_conformal_model2$pred[[1]]),
  ci_upper = unlist(nonamp_conformal_model2$up[[1]]),
  ci_lower = unlist(nonamp_conformal_model2$lo[[1]]),
  amp = Y_curves_shifted_mat[,156], 
  nonamp = resp_mat[,X_pred_nonamp_sim[2]]
) 

# Does the amputee exceed the bands?
ff_out_vec[2] <- with(non_amp_pred_band_data2_comp_plot, {
  sum(amp < ff_lower | amp > ff_upper)
})
ci_out_vec[2] <- with(non_amp_pred_band_data2_comp_plot, {
  sum(amp < ci_lower | amp > ci_upper)
})

# Comparison plot
conf_ff_comp2 <- ggplot() +
  geom_vline(xintercept = c(33 + 1/3, 66 + 2/3), color = "lightgray") +
  geom_ribbon(aes(x = 100*t_pts, ymin = ff_lower, ymax = ff_upper), 
              fill = "#56B4E9", alpha = 0.5, linewidth = 1.25,
              data = non_amp_pred_band_data2_comp_plot) +
  geom_ribbon(aes(x = 100*t_pts, ymin = ci_lower, ymax = ci_upper), 
              fill = NA, alpha = 0.5, color = "#D55E00",
              linewidth = 1.25, linetype = "dotdash",
              data = non_amp_pred_band_data2_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = nonamp), 
            color = "black", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data2_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = amp), 
            color = "black", linetype = "dashed", linewidth = 1.25,
            data = non_amp_pred_band_data2_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = ff_mean), 
            color = "#56B4E9", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data2_comp_plot) +
  scale_x_continuous(breaks = round(seq(0,100,(33 + 1/3)), 2)) +
  coord_cartesian(ylim = c(0, 1200), expand = FALSE) +
  theme_bw() +
  labs(x = "% Push-Off Phase", 
       y = "Force (N/kg)") + #, 
  # title = "Second Amputee") +
  theme( #legend.position = "none", 
    text = element_text(size = 20), 
    plot.margin = unit(c(.5, .9, .5, .5), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  annotate("label", 
           x = 35, y = 1100, 
           label = "Second Amputee", 
           size = 6) 

# ### Plot of the mean and band
# # Get the band
# non_amp_pred_band_data2 <- non_amp_pred_band2[[1]][[1]]
# # Mean
# par(mar = c(3,3,3,3))
# plot(grid, non_amp_pred_band_data2[,1],
#      type = "l", xlab = "T", ylab = "Front V Force",
#      xlim = c(0, .90),
#      ylim = c(0, max(Y_curves + 20)), 
#      main = paste0("95% Simultaneous Prediction Band for 2nd Amputee"))
# # Non-amputees
# for (i in nonamp_idx2) {
#   lines(x = grid,
#         y = resp_mat[, i],
#         col = "gray")
# }
# # Band
# lines(x = grid, y = non_amp_pred_band_data2[,1],
#       col = "black", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data2[,2],
#       col = "red", lty = "dashed", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data2[,3],
#       col = "red", lty = "dashed", lwd = 2)
# 
# # Amputee runner
# lines(x = grid, y = Y_curves_shifted_mat[,156], col = "green", lwd = 3)




## Third amputee runner
X_pred_amp3 <- X_pred_amp[3, c(2,3,4,6,10)]

# Non-amputee runners that are at least 23 years old and female
# The fifth runner of that^ list is the closest counterpart
nonamp_idx3 <- which(X_pred_nonamp$Sex == "male" &
                       X_pred_nonamp$Mass > 73 &
                       X_pred_nonamp$Mass < 76 &
                       X_pred_nonamp$TPush_V > 0.38 &
                       X_pred_nonamp$TPush_V < 0.40)[1]
# First amputee
new_dat_mat3 <- matrix(
  rep(c(X_pred_amp3$Mass[1], 
        X_pred_amp3$Age[1], 
        X_pred_amp3$Height[1], 
        1, 
        X_pred_amp3$TPush_V[1]), 
      each = n_sp), 
  nrow = n_sp
)
non_amp_pred_band3 <- predict_concurrent(
  concurrent_list = nonamp_model, 
  interval = "prediction", 
  new_dat = new_dat_mat3, 
  err_str = "t", 
  conf.level = 0.90,
  n_int = 3, 
  nu0_hat = "singh_df",
  mse_scalar = "ub"
)

### Plot of the mean and band
# Get the band
non_amp_pred_band_data3 <- non_amp_pred_band3[[1]][[1]]
non_amp_pred_band_data3_plot <- data.frame(
  t_pts = grid,
  mean = non_amp_pred_band_data3[,1],
  upper = non_amp_pred_band_data3[,2],
  lower = non_amp_pred_band_data3[,3],
  amp = Y_curves_shifted_mat[,157], 
  nonamp = resp_mat[,X_pred_nonamp_sim[3]]
) %>%
  pivot_longer(cols = -t_pts, 
               values_to = "Front_V", 
               names_to = "Estimate") 

amp_plot3 <- ggplot() +
  geom_line(aes(x = 100*t_pts, y = Front_V, 
                color = Estimate, linetype = Estimate), 
            data = non_amp_pred_band_data3_plot) +
  scale_color_manual(values = est_colors) +
  scale_linetype_manual(values = est_lines) +
  scale_x_continuous(breaks = seq(0,100,20)) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 20)) +
  labs(x = "% Push-Off Phase (realigned)", 
       y = "Front Vertical Force (N/kg)", 
       title = "Third Amputee") +
  coord_cartesian(ylim = c(0, 1250))
# title = paste0("95% Simultaneous Prediction Band to Predict ",
#                "3rd Amputee"),
# subtitle = paste0("Estimated for a ",
#                   "sprinter with Mass = 74.7 kg, Age = 25 ", 
#                   "years, Height = 1.91 m, Sex = Male, and ",
#                   "TPush_V = 0.388"))

### September 12th, 2023
## Creating a comparison plot of Conformal Inference vs Fast and Fair
non_amp_pred_band_data3_comp_plot <- data.frame(
  t_pts = grid,
  ff_mean = non_amp_pred_band_data3[,1],
  ff_upper = non_amp_pred_band_data3[,2],
  ff_lower = non_amp_pred_band_data3[,3],
  ci_mean = unlist(nonamp_conformal_model3$pred[[1]]),
  ci_upper = unlist(nonamp_conformal_model3$up[[1]]),
  ci_lower = unlist(nonamp_conformal_model3$lo[[1]]),
  amp = Y_curves_shifted_mat[,157], 
  nonamp = resp_mat[,X_pred_nonamp_sim[3]]
) 

# Does the amputee exceed the bands?
ff_out_vec[3] <- with(non_amp_pred_band_data3_comp_plot, {
  sum(amp < ff_lower | amp > ff_upper)
})
ci_out_vec[3] <- with(non_amp_pred_band_data3_comp_plot, {
  sum(amp < ci_lower | amp > ci_upper)
})

# Comparison plot
conf_ff_comp3 <- ggplot() +
  geom_vline(xintercept = c(33 + 1/3, 66 + 2/3), color = "lightgray") +
  geom_ribbon(aes(x = 100*t_pts, ymin = ff_lower, ymax = ff_upper), 
              fill = "#56B4E9", alpha = 0.5, linewidth = 1.25,
              data = non_amp_pred_band_data3_comp_plot) +
  geom_ribbon(aes(x = 100*t_pts, ymin = ci_lower, ymax = ci_upper), 
              fill = NA, alpha = 0.5, color = "#D55E00",
              linewidth = 1.25, linetype = "dotdash",
              data = non_amp_pred_band_data3_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = nonamp), 
            color = "black", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data3_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = amp), 
            color = "black", linetype = "dashed", linewidth = 1.25,
            data = non_amp_pred_band_data3_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = ff_mean), 
            color = "#56B4E9", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data3_comp_plot) +
  scale_x_continuous(breaks = round(seq(0,100,(33 + 1/3)), 2)) +
  coord_cartesian(ylim = c(0, 1200), expand = FALSE) +
  theme_bw() +
  labs(x = "% Push-Off Phase", 
       y = "Force (N/kg)") + #, 
  # title = "Third Amputee") +
  theme( #legend.position = "none", 
    text = element_text(size = 20), 
    plot.margin = unit(c(.5, .9, .5, .5), "cm"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  annotate("label", 
           x = 35, y = 1100, 
           label = "Third Amputee", 
           size = 6) 

# ### Plot of the mean and band
# # Get the band
# non_amp_pred_band_data3 <- non_amp_pred_band3[[1]][[1]]
# # Mean
# par(mar = c(3,3,3,3))
# plot(grid, non_amp_pred_band_data3[,1],
#      type = "l", xlab = "T", ylab = "Front V Force",
#      xlim = c(0, .90),
#      ylim = c(0, max(Y_curves + 20)), 
#      main = paste0("95% Simultaneous Prediction Band for 3rd Amputee"))
# # Non-amputees
# for (i in nonamp_idx3) {
#   lines(x = grid,
#         y = resp_mat[, i],
#         col = "gray")
# }
# # Band
# lines(x = grid, y = non_amp_pred_band_data3[,1],
#       col = "black", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data3[,2],
#       col = "red", lty = "dashed", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data3[,3],
#       col = "red", lty = "dashed", lwd = 2)
# 
# # Amputee runner
# lines(x = grid, y = Y_curves_shifted_mat[,157], col = "green", lwd = 3)




## Fourth amputee runner
X_pred_amp4 <- X_pred_amp[4, c(2,3,4,6,10)]

# Non-amputee runners that are at least 25 years old and female
# The fourth runner of that^ list is the closest counterpart
nonamp_idx4 <- which(X_pred_nonamp$Sex == "male" &
                       X_pred_nonamp$TPush_V > 0.43 &
                       X_pred_nonamp$Mass < 75)
# First amputee
new_dat_mat4 <- matrix(
  rep(c(X_pred_amp4$Mass[1], 
        X_pred_amp4$Age[1], 
        X_pred_amp4$Height[1], 
        1, 
        X_pred_amp4$TPush_V[1]), 
      each = n_sp), 
  nrow = n_sp
)
non_amp_pred_band4 <- predict_concurrent(
  concurrent_list = nonamp_model, 
  interval = "prediction", 
  new_dat = new_dat_mat4, 
  err_str = "t", 
  conf.level = 0.90,
  n_int = 3, 
  nu0_hat = "singh_df",
  mse_scalar = "ub"
)

### Plot of the mean and band
# Get the band
non_amp_pred_band_data4 <- non_amp_pred_band4[[1]][[1]]
non_amp_pred_band_data4_plot <- data.frame(
  t_pts = grid,
  mean = non_amp_pred_band_data4[,1],
  upper = non_amp_pred_band_data4[,2],
  lower = non_amp_pred_band_data4[,3],
  amp = Y_curves_shifted_mat[,158], 
  nonamp = resp_mat[,X_pred_nonamp_sim[4]]
) %>%
  pivot_longer(cols = -t_pts, 
               values_to = "Front_V", 
               names_to = "Estimate") 

amp_plot4 <- ggplot() +
  geom_line(aes(x = 100*t_pts, y = Front_V, 
                color = Estimate, linetype = Estimate), 
            data = non_amp_pred_band_data4_plot) +
  scale_color_manual(values = est_colors) +
  scale_linetype_manual(values = est_lines) +
  scale_x_continuous(breaks = seq(0,100,20)) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 20)) +
  labs(x = "% Push-Off Phase (realigned)", 
       y = "Front Vertical Force (N/kg)", 
       title = "Fourth Amputee") +
  coord_cartesian(ylim = c(0, 1250))
# title = paste0("95% Simultaneous Prediction Band to Predict ",
#                "4th Amputee"),
# subtitle = paste0("Estimated for a ",
#                   "sprinter with Mass = 69.7 kg, Age = 27 ", 
#                   "years, Height = 1.87 m, Sex = Male, and ",
#                   "TPush_V = 0.466"))

### September 12th, 2023
## Creating a comparison plot of Conformal Inference vs Fast and Fair
non_amp_pred_band_data4_comp_plot <- data.frame(
  t_pts = grid,
  ff_mean = non_amp_pred_band_data4[,1],
  ff_upper = non_amp_pred_band_data4[,2],
  ff_lower = non_amp_pred_band_data4[,3],
  ci_mean = unlist(nonamp_conformal_model4$pred[[1]]),
  ci_upper = unlist(nonamp_conformal_model4$up[[1]]),
  ci_lower = unlist(nonamp_conformal_model4$lo[[1]]),
  amp = Y_curves_shifted_mat[,158], 
  nonamp = resp_mat[,X_pred_nonamp_sim[4]]
) 

# Does the amputee exceed the bands?
ff_out_vec[4] <- with(non_amp_pred_band_data4_comp_plot, {
  sum(amp < ff_lower | amp > ff_upper)
})
ci_out_vec[4] <- with(non_amp_pred_band_data4_comp_plot, {
  sum(amp < ci_lower | amp > ci_upper)
})

# Comparison plot
conf_ff_comp4 <- ggplot() +
  geom_vline(xintercept = c(33 + 1/3, 66 + 2/3), color = "lightgray") +
  geom_ribbon(aes(x = 100*t_pts, ymin = ff_lower, ymax = ff_upper), 
              fill = "#56B4E9", alpha = 0.5, linewidth = 1.25,
              data = non_amp_pred_band_data4_comp_plot) +
  geom_ribbon(aes(x = 100*t_pts, ymin = ci_lower, ymax = ci_upper), 
              fill = NA, alpha = 0.5, color = "#D55E00",
              linewidth = 1.25, linetype = "dotdash",
              data = non_amp_pred_band_data4_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = nonamp), 
            color = "black", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data4_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = amp), 
            color = "black", linetype = "dashed", linewidth = 1.25,
            data = non_amp_pred_band_data4_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = ff_mean), 
            color = "#56B4E9", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data4_comp_plot) +
  scale_x_continuous(breaks = round(seq(0,100,(33 + 1/3)), 2)) +
  coord_cartesian(ylim = c(0, 1200), expand = FALSE) +
  theme_bw() +
  labs(x = "% Push-Off Phase", 
       y = "Force (N/kg)") + #, 
  # title = "Fourth Amputee") +
  theme( #legend.position = "none", 
    text = element_text(size = 20), 
    plot.margin = unit(c(.5, .9, .5, .5), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  annotate("label", 
           x = 35, y = 1100, 
           label = "Fourth Amputee", 
           size = 6) 

# ### Plot of the mean and band
# # Get the band
# non_amp_pred_band_data4 <- non_amp_pred_band4[[1]][[1]]
# # Mean
# par(mar = c(3,3,3,3))
# plot(grid, non_amp_pred_band_data4[,1],
#      type = "l", xlab = "T", ylab = "Front V Force",
#      xlim = c(0, .90),
#      ylim = c(0, max(Y_curves + 20)), 
#      main = paste0("95% Simultaneous Prediction Band for 4th Amputee"))
# # Non-amputees
# for (i in nonamp_idx4) {
#   lines(x = grid,
#         y = resp_mat[, i],
#         col = "gray")
# }
# # Band
# lines(x = grid, y = non_amp_pred_band_data4[,1],
#       col = "black", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data4[,2],
#       col = "red", lty = "dashed", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data4[,3],
#       col = "red", lty = "dashed", lwd = 2)
# 
# # Amputee runner
# lines(x = grid, y = Y_curves_shifted_mat[,158], col = "green", lwd = 3)



## Fifth amputee runner
X_pred_amp5 <- X_pred_amp[5, c(2,3,4,6,10)]

# Non-amputee runners that are at least 25 years old and female
# The fourth runner of that^ list is the closest counterpart
nonamp_idx5 <- which(X_pred_nonamp$Mass > 80 &
                       X_pred_nonamp$Mass < 82 &
                       X_pred_nonamp$Sex == "male" &
                       X_pred_nonamp$TPush_V > 0.375)
# Fifth amputee
new_dat_mat5 <- matrix(
  rep(c(X_pred_amp5$Mass[1], 
        X_pred_amp5$Age[1], 
        X_pred_amp5$Height[1], 
        1, 
        X_pred_amp5$TPush_V[1]), 
      each = n_sp), 
  nrow = n_sp
)
non_amp_pred_band5 <- predict_concurrent(
  concurrent_list = nonamp_model, 
  interval = "prediction", 
  new_dat = new_dat_mat5, 
  err_str = "t", 
  conf.level = 0.90,
  n_int = 3, 
  nu0_hat = "singh_df",
  mse_scalar = "ub"
)

### Plot of the mean and band
# Get the band
non_amp_pred_band_data5 <- non_amp_pred_band5[[1]][[1]]
non_amp_pred_band_data5_plot <- data.frame(
  t_pts = grid,
  mean = non_amp_pred_band_data5[,1],
  upper = non_amp_pred_band_data5[,2],
  lower = non_amp_pred_band_data5[,3],
  amp = Y_curves_shifted_mat[,159], 
  nonamp = resp_mat[,X_pred_nonamp_sim[5]]
) %>%
  pivot_longer(cols = -t_pts, 
               values_to = "Front_V", 
               names_to = "Estimate") 

amp_plot5 <- ggplot() +
  geom_line(aes(x = 100*t_pts, y = Front_V, 
                color = Estimate, linetype = Estimate), 
            data = non_amp_pred_band_data5_plot) +
  scale_color_manual(values = est_colors) +
  scale_linetype_manual(values = est_lines) +
  scale_x_continuous(breaks = seq(0,100,20)) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 20)) +
  labs(x = "% Push-Off Phase (realigned)", 
       y = "Front Vertical Force (N/kg)", 
       title = "Fifth Amputee") +
  coord_cartesian(ylim = c(0, 1250))
# title = paste0("95% Simultaneous Prediction Band to Predict ",
#                "5th Amputee"),
# subtitle = paste0("Estimated for a ",
#                   "sprinter with Mass = 80.2 kg, Age = 30 ", 
#                   "years, Height = 1.81 m, Sex = Male, and ",
#                   "TPush_V = 0.441"))

### September 12th, 2023
## Creating a comparison plot of Conformal Inference vs Fast and Fair
non_amp_pred_band_data5_comp_plot <- data.frame(
  t_pts = grid,
  ff_mean = non_amp_pred_band_data5[,1],
  ff_upper = non_amp_pred_band_data5[,2],
  ff_lower = non_amp_pred_band_data5[,3],
  ci_mean = unlist(nonamp_conformal_model5$pred[[1]]),
  ci_upper = unlist(nonamp_conformal_model5$up[[1]]),
  ci_lower = unlist(nonamp_conformal_model5$lo[[1]]),
  amp = Y_curves_shifted_mat[,159], 
  nonamp = resp_mat[,X_pred_nonamp_sim[5]]
) 

# Does the amputee exceed the bands?
ff_out_vec[5] <- with(non_amp_pred_band_data5_comp_plot, {
  sum(amp < ff_lower | amp > ff_upper)
})
ci_out_vec[5] <- with(non_amp_pred_band_data5_comp_plot, {
  sum(amp < ci_lower | amp > ci_upper)
})

# Comparison plot
conf_ff_comp5 <- ggplot() +
  geom_vline(xintercept = c(33 + 1/3, 66 + 2/3), color = "lightgray") +
  geom_ribbon(aes(x = 100*t_pts, ymin = ff_lower, ymax = ff_upper), 
              fill = "#56B4E9", alpha = 0.5, linewidth = 1.25,
              data = non_amp_pred_band_data5_comp_plot) +
  geom_ribbon(aes(x = 100*t_pts, ymin = ci_lower, ymax = ci_upper), 
              fill = NA, alpha = 0.5, color = "#D55E00",
              linewidth = 1.25, linetype = "dotdash",
              data = non_amp_pred_band_data5_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = nonamp), 
            color = "black", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data5_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = amp), 
            color = "black", linetype = "dashed", linewidth = 1.25,
            data = non_amp_pred_band_data5_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = ff_mean), 
            color = "#56B4E9", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data5_comp_plot) +
  scale_x_continuous(breaks = round(seq(0,100,(33 + 1/3)), 2)) +
  coord_cartesian(ylim = c(0, 1200), expand = FALSE) +
  theme_bw() +
  labs(x = "% Push-Off Phase", 
       y = "Force (N/kg)") + #, 
  # title = "Fifth Amputee") +
  theme( #legend.position = "none", 
    text = element_text(size = 20), 
    plot.margin = unit(c(.5, .9, .5, .5), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  annotate("label", 
           x = 35, y = 1100, 
           label = "Fifth Amputee", 
           size = 6)

# ### Plot of the mean and band
# # Get the band
# non_amp_pred_band_data5 <- non_amp_pred_band5[[1]][[1]]
# # Mean
# par(mar = c(3,3,3,3))
# plot(grid, non_amp_pred_band_data5[,1],
#      type = "l", xlab = "T", ylab = "Front V Force",
#      xlim = c(0, .90),
#      ylim = c(0, max(Y_curves + 20)), 
#      main = paste0("95% Simultaneous Prediction Band for 5th Amputee"))
# # Non-amputees
# for (i in nonamp_idx5) {
#   lines(x = grid,
#         y = resp_mat[, i],
#         col = "gray")
# }
# # Band
# lines(x = grid, y = non_amp_pred_band_data5[,1],
#       col = "black", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data5[,2],
#       col = "red", lty = "dashed", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data5[,3],
#       col = "red", lty = "dashed", lwd = 2)
# 
# # Amputee runner
# lines(x = grid, y = Y_curves_shifted_mat[,159], col = "green", lwd = 3)



## Sixth amputee runner
X_pred_amp6 <- X_pred_amp[6, c(2,3,4,6,10)]

# Non-amputee runners that are at least 25 years old and female
# The fourth runner of that^ list is the closest counterpart
nonamp_idx6 <- which(X_pred_nonamp$Age > 23 &
                       X_pred_nonamp$Age < 27 &
                       X_pred_nonamp$Sex == "male" &
                       X_pred_nonamp$Mass > 87 &
                       X_pred_nonamp$Mass < 92)
# First amputee
new_dat_mat6 <- matrix(
  rep(c(X_pred_amp6$Mass[1], 
        X_pred_amp6$Age[1], 
        X_pred_amp6$Height[1], 
        1, 
        X_pred_amp6$TPush_V[1]), 
      each = n_sp), 
  nrow = n_sp
)
non_amp_pred_band6 <- predict_concurrent(
  concurrent_list = nonamp_model, 
  interval = "prediction", 
  new_dat = new_dat_mat6, 
  err_str = "t", 
  conf.level = 0.90,
  n_int = 3, 
  nu0_hat = "singh_df",
  mse_scalar = "ub"
)

### Plot of the mean and band
# Get the band
non_amp_pred_band_data6 <- non_amp_pred_band6[[1]][[1]]
non_amp_pred_band_data6_plot <- data.frame(
  t_pts = grid,
  mean = non_amp_pred_band_data6[,1],
  upper = non_amp_pred_band_data6[,2],
  lower = non_amp_pred_band_data6[,3],
  amp = Y_curves_shifted_mat[,160], 
  nonamp = resp_mat[,X_pred_nonamp_sim[6]]
) %>%
  pivot_longer(cols = -t_pts, 
               values_to = "Front_V", 
               names_to = "Estimate") 

amp_plot6 <- ggplot() +
  geom_line(aes(x = 100*t_pts, y = Front_V, 
                color = Estimate, linetype = Estimate), 
            data = non_amp_pred_band_data6_plot) +
  scale_color_manual(values = est_colors) +
  scale_linetype_manual(values = est_lines) +
  scale_x_continuous(breaks = seq(0,100,20)) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 20)) +
  labs(x = "% Push-Off Phase (realigned)", 
       y = "Front Vertical Force (N/kg)", 
       title = "Sixth Amputee") +
  coord_cartesian(ylim = c(0, 1250))
# title = paste0("95% Simultaneous Prediction Band to Predict ",
#                "6th Amputee"),
# subtitle = paste0("Estimated for a ",
#                   "sprinter with Mass = 89.1 kg, Age = 25 ", 
#                   "years, Height = 1.97 m, Sex = Male, and ",
#                   "TPush_V = 0.373"))

### September 12th, 2023
## Creating a comparison plot of Conformal Inference vs Fast and Fair
non_amp_pred_band_data6_comp_plot <- data.frame(
  t_pts = grid,
  ff_mean = non_amp_pred_band_data6[,1],
  ff_upper = non_amp_pred_band_data6[,2],
  ff_lower = non_amp_pred_band_data6[,3],
  ci_mean = unlist(nonamp_conformal_model6$pred[[1]]),
  ci_upper = unlist(nonamp_conformal_model6$up[[1]]),
  ci_lower = unlist(nonamp_conformal_model6$lo[[1]]),
  amp = Y_curves_shifted_mat[,160], 
  nonamp = resp_mat[,X_pred_nonamp_sim[6]]
) 

# Does the amputee exceed the bands?
ff_out_vec[6] <- with(non_amp_pred_band_data6_comp_plot, {
  sum(amp < ff_lower | amp > ff_upper)
})
ci_out_vec[6] <- with(non_amp_pred_band_data6_comp_plot, {
  sum(amp < ci_lower | amp > ci_upper)
})

# Comparison plot
conf_ff_comp6 <- ggplot() +
  geom_vline(xintercept = c(33 + 1/3, 66 + 2/3), color = "lightgray") +
  geom_ribbon(aes(x = 100*t_pts, ymin = ff_lower, ymax = ff_upper), 
              fill = "#56B4E9", alpha = 0.5, linewidth = 1.25,
              data = non_amp_pred_band_data6_comp_plot) +
  geom_ribbon(aes(x = 100*t_pts, ymin = ci_lower, ymax = ci_upper), 
              fill = NA, alpha = 0.5, color = "#D55E00",
              linewidth = 1.25, linetype = "dotdash",
              data = non_amp_pred_band_data6_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = nonamp), 
            color = "black", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data6_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = amp), 
            color = "black", linetype = "dashed", linewidth = 1.25,
            data = non_amp_pred_band_data6_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = ff_mean), 
            color = "#56B4E9", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data6_comp_plot) +
  scale_x_continuous(breaks = round(seq(0,100,(33 + 1/3)), 2)) +
  coord_cartesian(ylim = c(0, 1200), expand = FALSE) +
  theme_bw() +
  labs(x = "% Push-Off Phase", 
       y = "Force (N/kg)") + #, 
  # title = "Sixth Amputee") +
  theme( #legend.position = "none", 
    text = element_text(size = 20), 
    plot.margin = unit(c(.5, .9, .5, .5), "cm"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  annotate("label", 
           x = 35, y = 1100, 
           label = "Sixth Amputee", 
           size = 6) 

# ### Plot of the mean and band
# # Get the band
# non_amp_pred_band_data6 <- non_amp_pred_band6[[1]][[1]]
# # Mean
# par(mar = c(3,3,3,3))
# plot(grid, non_amp_pred_band_data6[,1],
#      type = "l", xlab = "T", ylab = "Front V Force",
#      xlim = c(0, .90),
#      ylim = c(0, max(Y_curves + 20)), 
#      main = paste0("95% Simultaneous Prediction Band for 6th Amputee"))
# # Non-amputees
# for (i in nonamp_idx6) {
#   lines(x = grid,
#         y = resp_mat[, i],
#         col = "gray")
# }
# # Band
# lines(x = grid, y = non_amp_pred_band_data6[,1],
#       col = "black", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data6[,2],
#       col = "red", lty = "dashed", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data6[,3],
#       col = "red", lty = "dashed", lwd = 2)
# 
# # Amputee runner
# lines(x = grid, y = Y_curves_shifted_mat[,160], col = "green", lwd = 3)



## Seventh amputee runner
X_pred_amp7 <- X_pred_amp[7, c(2,3,4,6,10)]

# Non-amputee runners that are at least 25 years old and female
# The fourth runner of that^ list is the closest counterpart
nonamp_idx7 <- which(X_pred_nonamp$Sex == "male" &
                       (X_pred_nonamp$Mass < 72 &
                          X_pred_nonamp$Mass > 70) &
                       X_pred_nonamp$TPush_V > 0.391)
# First amputee
new_dat_mat7 <- matrix(
  rep(c(X_pred_amp7$Mass[1], 
        X_pred_amp7$Age[1], 
        X_pred_amp7$Height[1], 
        1, 
        X_pred_amp7$TPush_V[1]), 
      each = n_sp), 
  nrow = n_sp
)
non_amp_pred_band7 <- predict_concurrent(
  concurrent_list = nonamp_model, 
  interval = "prediction", 
  new_dat = new_dat_mat7, 
  err_str = "t", 
  conf.level = 0.90,
  n_int = 3, 
  nu0_hat = "singh_df",
  mse_scalar = "ub"
)

### Plot of the mean and band
# Get the band
non_amp_pred_band_data7 <- non_amp_pred_band7[[1]][[1]]
non_amp_pred_band_data7_plot <- data.frame(
  t_pts = grid,
  mean = non_amp_pred_band_data7[,1],
  upper = non_amp_pred_band_data7[,2],
  lower = non_amp_pred_band_data7[,3],
  amp = Y_curves_shifted_mat[,161], 
  nonamp = resp_mat[,X_pred_nonamp_sim[7]]
) %>%
  pivot_longer(cols = -t_pts, 
               values_to = "Front_V", 
               names_to = "Estimate") 

amp_plot7 <- ggplot() +
  geom_line(aes(x = 100*t_pts, y = Front_V, 
                color = Estimate, linetype = Estimate), 
            data = non_amp_pred_band_data7_plot) +
  geom_vline(xintercept = c(33 + 1/3, 66 + 2/3)) +
  scale_color_manual(values = est_colors) +
  scale_linetype_manual(values = est_lines) +
  scale_x_continuous(breaks = seq(0,100,(33 + 1/3))) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 20)) +
  labs(x = "% Push-Off Phase (realigned)", 
       y = "Force (N/kg)", 
       title = "Seventh Amputee") +
  coord_cartesian(ylim = c(0, 1250))
# title = paste0("95% Simultaneous Prediction Band to Predict ",
#                "7th Amputee"),
# subtitle = paste0("Estimated for a ",
#                   "sprinter with Mass = 71.0 kg, Age = 31 ", 
#                   "years, Height = 1.78 m, Sex = Male, and ",
#                   "TPush_V = 0.416"))


### September 12th, 2023
## Creating a comparison plot of Conformal Inference vs Fast and Fair
non_amp_pred_band_data7_comp_plot <- data.frame(
  t_pts = grid,
  ff_mean = non_amp_pred_band_data7[,1],
  ff_upper = non_amp_pred_band_data7[,2],
  ff_lower = non_amp_pred_band_data7[,3],
  ci_mean = unlist(nonamp_conformal_model7$pred[[1]]),
  ci_upper = unlist(nonamp_conformal_model7$up[[1]]),
  ci_lower = unlist(nonamp_conformal_model7$lo[[1]]),
  amp = Y_curves_shifted_mat[,161], 
  nonamp = resp_mat[,X_pred_nonamp_sim[7]]
) 

# Does the amputee exceed the bands?
ff_out_vec[7] <- with(non_amp_pred_band_data7_comp_plot, {
  sum(amp < ff_lower | amp > ff_upper)
})
ci_out_vec[7] <- with(non_amp_pred_band_data7_comp_plot, {
  sum(amp < ci_lower | amp > ci_upper)
})

# Elements for manual legend
colors <- c("Fast n Fair" = "#56B4E9", 
            "Conformal Inference" = "#D55E00", 
            "Seventh Amputee" = "black")
linetypes <- c("Fast n Fair" = "solid", 
               "Conformal Inference" = "dotdash", 
               "Seventh Amputee" = "dashed")

## Recreates Figure 1 in Creutzinger, Liebl, and Sharp (2024+)
# Comparison plot
conf_ff_comp7_fig1 <- ggplot() +
  geom_vline(xintercept = c(33 + 1/3, 66 + 2/3), color = "lightgray") +
  geom_ribbon(aes(x = 100*t_pts, ymin = ff_lower, ymax = ff_upper), 
              fill = "#56B4E9", alpha = 0.5, linewidth = 1.25,
              data = non_amp_pred_band_data7_comp_plot, 
              key_glyph = "blank") +
  geom_ribbon(aes(x = 100*t_pts, ymin = ci_lower, ymax = ci_upper, 
                  color = "Conformal Inference",
                  linetype = "Conformal Inference"),
              fill = NA, alpha = 0.5,  
              linewidth = 1.25, 
              data = non_amp_pred_band_data7_comp_plot, 
              key_glyph = "blank") +
  # geom_line(aes(x = 100*t_pts, y = nonamp), 
  #               color = "black", linetype = "solid", linewidth = 1.25,
  #           data = non_amp_pred_band_data7_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = amp, 
                color = "Seventh Amputee", 
                linetype = "Seventh Amputee"), 
            linewidth = 1.25,
            data = non_amp_pred_band_data7_comp_plot, 
            key_glyph = "blank") +
  geom_line(aes(x = 100*t_pts, y = ff_mean, 
                color = "Fast n Fair Prediction Band", 
                linetype = "Fast n Fair Prediction Band"), 
            linewidth = 1.25,
            data = non_amp_pred_band_data7_comp_plot) +
  scale_x_continuous(breaks = round(seq(0,100,(33 + 1/3)), 2)) +
  coord_cartesian(ylim = c(0, 1200), expand = FALSE) +
  theme_bw() +
  scale_color_manual(values = colors) + 
  scale_linetype_manual(values = linetypes) +
  labs(x = "% Push-Off Phase", 
       y = "Force (N/kg)", 
       color = "", 
       linetype = "") +  
  # title = "Seventh Amputee") +
  theme(text = element_text(size = 20), 
        plot.margin = unit(c(1, 1, 1, 1), "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.28, 0.9), 
        legend.title = element_blank(), 
        legend.box.background = element_rect(color = "black"), 
        legend.spacing.y = unit(0, "mm"), 
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 16))
# annotate("label", 
#          x = 35, y = 1100, 
#          label = "Seventh Amputee", 
#          size = 6)

conf_ff_comp7 <- ggplot() +
  geom_vline(xintercept = c(33 + 1/3, 66 + 2/3), color = "lightgray") +
  geom_ribbon(aes(x = 100*t_pts, ymin = ff_lower, ymax = ff_upper), 
              fill = "#56B4E9", alpha = 0.5, linewidth = 1.25,
              data = non_amp_pred_band_data7_comp_plot) +
  geom_ribbon(aes(x = 100*t_pts, ymin = ci_lower, ymax = ci_upper), 
              fill = NA, alpha = 0.5, color = "#D55E00",
              linewidth = 1.25, linetype = "dotdash",
              data = non_amp_pred_band_data7_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = nonamp), 
            color = "black", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data7_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = amp), 
            color = "black", linetype = "dashed", linewidth = 1.25,
            data = non_amp_pred_band_data7_comp_plot) +
  geom_line(aes(x = 100*t_pts, y = ff_mean), 
            color = "#56B4E9", linetype = "solid", linewidth = 1.25,
            data = non_amp_pred_band_data7_comp_plot) +
  scale_x_continuous(breaks = round(seq(0,100,(33 + 1/3)), 2)) +
  coord_cartesian(ylim = c(0, 1200), expand = FALSE) +
  theme_bw() +
  labs(x = "% Push-Off Phase", 
       y = "Force (N/kg)") + #, 
  # title = "Sixth Amputee") +
  theme( #legend.position = "none", 
    text = element_text(size = 20), 
    plot.margin = unit(c(.5, .9, .5, .5), "cm"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  annotate("label", 
           x = 35, y = 1100, 
           label = "Seventh Amputee", 
           size = 6) 


# ### Plot of the mean and band
# # Get the band
# non_amp_pred_band_data7 <- non_amp_pred_band7[[1]][[1]]
# # Mean
# par(mar = c(3,3,3,3))
# plot(grid, non_amp_pred_band_data7[,1],
#      type = "l", xlab = "T", ylab = "Front V Force",
#      xlim = c(0, .90),
#      ylim = c(0, max(Y_curves + 20)), 
#      main = paste0("95% Simultaneous Prediction Band for 7th Amputee"))
# # Non-amputees
# for (i in nonamp_idx7) {
#   lines(x = grid,
#         y = resp_mat[, i],
#         col = "gray")
# }
# # Band
# lines(x = grid, y = non_amp_pred_band_data7[,1],
#       col = "black", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data7[,2],
#       col = "red", lty = "dashed", lwd = 2)
# lines(x = grid, y = non_amp_pred_band_data7[,3],
#       col = "red", lty = "dashed", lwd = 2)
# 
# # Amputee runner
# lines(x = grid, y = Y_curves_shifted_mat[,161], col = "green", lwd = 3)

## Arrange all the plots into one figure
# Make a "figure" of demographics table
X_mod_amp <- X_pred_amp[, c(2,3,4,6,10)] %>%
  rename("Push-Time" = "TPush_V")

demo_tab <- tableGrob(X_mod_amp, 
                      theme = ttheme_default(base_size = 16))
demo_fig <- ggplot() + 
  geom_point(x = 1:10, y = 1:10) +
  annotation_custom(tableGrob(X_mod_amp,
                              theme = ttheme_default(base_size = 16)),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  labs(title = "Amputee Demographics \n") +
  theme_bw() +
  theme(
    title = element_text(hjust = 1, vjust = -0.5),
    text = element_text(size = 20),
    panel.spacing = unit(1, "cm"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0.5, 0.9, 1, 2), "cm")
  )

grid.arrange(amp_plot1, amp_plot2, amp_plot3, amp_plot4, amp_plot5, 
             amp_plot6, amp_plot7, demo_fig, ncol = 2)

grid.arrange(conf_ff_comp1, conf_ff_comp2, conf_ff_comp3, 
             conf_ff_comp4, conf_ff_comp5, conf_ff_comp6, 
             conf_ff_comp7, demo_fig, ncol = 2)
#####################################################################