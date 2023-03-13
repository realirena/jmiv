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

##---------------------------------------
# import the data 
##---------------------------------------
## we have already excluded all of the people who were on hormone therapy i
hormone_data <- data.frame(read.csv("/home/irena/gsra/swan_sims/data_analysis/visceral_fat/hormone_whr_10132022.csv"))
bm_data <- read.csv("/home/irena/gsra/swan_sims/data_analysis/visceral_fat/whr_data_10132022.csv")

seed = 607

Q = 2 # no of hormones
P = 2 # no of basis functions 
I = length(unique(bm_data$new_id))

hormone_mat <- data.matrix(hormone_data[, c("e2_loess_resid","fsh_loess_resid")])

cov_mat <- bm_data[,c("first_waist","black_ind", "chin_ind" ,"jpn_ind", "his_ind", "sports_cat2", "sports_cat3", "sports_cat4", "sports_cat5")]
##---------------------------------------:
## compile and run stan model: 
##---------------------------------------

joint_model <- stan_model("/home/irena/gsra/swan_sims/data_analysis/visceral_fat/joint_model_mixture_scale_w_cor.stan")

sim_out <- sampling(joint_model, sample_file=paste0('/home/irena/gsra/swan_sims/data_analysis/visceral_fat/1013_waist_model_samples.csv'),
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
                                                         bm_outcome = bm_data$waist_rate_change,
                                                         simulate = 0
                                             ),
                         seed=seed
                        )
