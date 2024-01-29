#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### August 16th, 2023                                             ###
###                                                               ###
###   Simulation study results.                                   ###
#####################################################################



## Starting steps
#####################################################################
# Install ffscb from forked package
require(devtools)
install_github("creutzml/ffscb")

# Packages
library(mvtnorm)
library(tidyverse)
library(fda)
library(ffscb)
library(conformalInference.fd)
#####################################################################



### Analyze the results:
#####################################################################
## Data Directory
dir_path <- file.path(here::here(), "Simulations/data")

## Load in the simulation results
# Stationary Matern covariance results
load(file.path(dir_path, "concurrent_sim_st.RData"))
# Non-stationary Matern covariance results
# load(file.path(dir_path, "concurrent_sim_nonst.RData"))

# Recreates results for Tables 1 and 3 (if stationary results are 
# loaded), and Tables 2 and 4 (if non-stationary results are loaded)
## Modify variables appropriately and create a summary table:
sim_results_df2 <- sim_results_df %>%
  mutate(out = as.logical(out),
         out1 = as.logical(out1),
         out2 = as.logical(out2),
         out3 = as.logical(out3),
         # out4 = as.logical(out4),
         # out5 = as.logical(out5),
         band_score = as.numeric(band_score), 
         max_band_width = as.numeric(max_band_width)) %>%
  group_by(n, dfs, method, prediction) %>%
  summarise(covg = 1 - mean(out),
            covg_int1 = 1 - mean(out1),
            covg_int2 = 1 - mean(out2),
            covg_int3 = 1 - mean(out3),
            # covg_int4 = 1 - mean(out4),
            # covg_int5 = 1 - mean(out5),
            mean_band_score = mean(band_score, na.rm = T), 
            mean_max_band_width = mean(max_band_width, 
                                       na.rm = T)) %>%
  mutate(covg = format(round(covg, 2), 
                       nsmall = 2),
         covg_int1 = format(round(covg_int1, 2), 
                       nsmall = 2),
         covg_int2 = format(round(covg_int2, 2), 
                       nsmall = 2),
         covg_int3 = format(round(covg_int3, 2), 
                       nsmall = 2),
         # covg_int4 = format(round(covg_int4, 2), 
                       # nsmall = 2),
         # covg_int5 = format(round(covg_int5, 2), 
                       # nsmall = 2),
         mean_band_score = format(round(mean_band_score, 2), 
                                  nsmall = 2),
         mean_max_band_width = format(round(mean_max_band_width, 2), 
                                      nsmall = 2))
  # pivot_wider(names_from = c(method, prediction), 
  #             values_from = c(covg:mean_band_width))

#####################################################################