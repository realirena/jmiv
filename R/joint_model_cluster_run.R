##---------------------------------------
## libraries:
##---------------------------------------
rm(list=ls())
library(rstan)
library(dplyr)
library(ggplot2)
library(reshape2)
library(splines)
library(MASS)
library(StanHeaders)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## set seed 
seed = 607
##---------------------------------------
# import the data 
##---------------------------------------

## predictor data (longitudinal)
hormone_data <- data.frame(read.csv("hormone_whr_10132022.csv"))

## outcome data (cross-sectional)
bm_data <- read.csv("whr_data_10132022.csv")

Q = 2 # no of hormones
P = 2 # no of basis functions for each marker 
I = length(unique(bm_data$new_id)) ## number of individuals 

## matrix for the predictor (markers)
hormone_mat <- data.matrix(hormone_data[, c("e2_loess_resid","fsh_loess_resid")])

## covariates in the outcome 
cov_mat <- bm_data[,c("first_waist","black_ind", "chin_ind" ,"jpn_ind", "his_ind", "sports_cat2", "sports_cat3", "sports_cat4", "sports_cat5")]
##---------------------------------------:
## compile and run stan model: 
##---------------------------------------

joint_model <- stan_model("joint_model_mixture_scale_w_cov.stan")

sim_out <- sampling(joint_model, 
                         iter = 2000,
                         warmup=1000, #BURN IN
                         chains = 4,
                         control = list(max_treedepth = 30,
                                        adapt_delta=0.95),
                                             data = list(Q = Q,
                                                         P = P,
                                                         I = I,
                                                         N = nrow(hormone_mat),
                                                         id=hormone_data$new_id,
                                                         f_time_fmp = cbind(1, hormone_data$time_fmp),
                                                         hormones =hormone_mat,
                                                         n_cov = 9,
                                                         K = 2,
                                                         dirichlet_alpha=1,
                                                         other_covariates = cov_mat,
                                                         outcome = bm_data$waist_rate_change,
                                                         simulate = 0),
                    seed=seed
                        )
