## cluster file 

rm(list=ls())
library(rstan)
library(dplyr)
library(reshape2)
library(splines)
library(MASS)
options(mc.cores = parallel::detectCores(logical= FALSE))
# rstan_options(auto_write = TRUE)
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

taskid <- as.numeric(slurm_arrayid)

#setwd("/home/irena/")
## set up the data simulation parameters:
seed = 511+taskid*100
set.seed(seed)

Q = 2 # no of hormones
P = 2 # no of basis functions 
I = 1000 # no of subjects ## increase subjects 

S <- lapply(1:I, function(x){
  q <- c(rnorm(1, 0, 0.75)/2, rnorm(1, 0.5, 0.5)/2)
  Q_mat <- diag(exp(q))
  R_mat = diag(Q)
  t = rbeta(1,1,5)
  r = 2*t-1
  R_mat[1,Q] <- r
  R_mat[Q,1] <- R_mat[1,Q]
  return(Q_mat %*% R_mat %*%Q_mat)
})

## import estimated betas from the hormone data: 
beta <- list()
beta[[1]] <- c(0,2)
beta[[2]] <- c(2,1)

## import estimated Sigmas from the hormone data: 
Sigma <- list()
Sigma[[1]] <- matrix(c(1, -0.05, -0.05, 1), ncol=P, nrow=P)
Sigma[[2]] <-  matrix(c(1, -0.1, -0.1, 0.5), ncol=P, nrow=P)

#---------------------------------------
### simulate the time to FMP and the hormones now: 
##---------------------------------------
Ti = sample(6:15, I, replace=T) ##no visits per subject

## placeholder for the simulated time to FMPs
time_fmp = rep(0, sum(Ti)) 

ids <- unlist(lapply(seq_along(Ti), function(i) rep(i, Ti[i])))

## simulate the time variable
for(i in 1:I){
  if(i==1){
    time_fmp = (seq(1:Ti[i]) + (sample(c(-9:-1),1)) + rnorm(Ti[i],0,0.5)) ## addnoise
    # apply basis function to time fmp 
    #f_time_fmp = bs(time_fmp, df=P)
    f_time_fmp = cbind(1, time_fmp)
  } else {
    time_fmp_i = seq(1:Ti[i]) + (sample(c(-9:-1),1)) + rnorm(Ti[i],0,0.5)
    time_fmp = c(time_fmp, time_fmp_i)
    #f_time_fmp = rbind(f_time_fmp, bs(time_fmp_i, df=P))
    f_time_fmp = rbind(f_time_fmp, cbind(1, time_fmp_i))
  }
}

### Generate the individual level coefficients 
B  <- lapply(1:I, function(x){
  mat = matrix(ncol=Q, nrow=P)
  for(q in 1:Q){
    mat[,q] <- mvrnorm(1, beta[[q]], Sigma[[q]])
  }
  return(mat)
})

## compute the means of the hormones
sim_mu <- t(sapply(seq(nrow(f_time_fmp)), function(i){
  B_i <- B[[ids[i]]]
  return(f_time_fmp[i,] %*% B_i)
}))

### generate the hormone data 
# sim_hormones <- t(sapply(seq_along(ids), function(i){mvrnorm(1, sim_mu[i,], S[[ids[i]]])}))
sim_hormones <- t(sapply(seq_along(ids), function(i){mvrnorm(1, sim_mu[i,], S[[ids[i]]])}))


## unlist B into a matrix with nrow = I and ncol = P*Q
B_design <- matrix(unlist(B), ncol=P*Q, byrow=T)

## make each S_i an upper triangular matrix 
S_lower <- lapply(1:I, function(i){
  mat = S[[i]]
  mat[upper.tri(mat)==T] <- 0
  return(mat)
})
## unroll S into a design matrix: 
S_design <- matrix(unlist(S_lower), ncol=Q*Q, byrow=T)


## set the true coefficient parameters for the means and variances 
alpha <- c(-3, -3, -3, 3)
gamma <- c(2,-1, 0,-2)

## get the mean of the outcome variable: 
mu_outcome = (B_design%*%alpha + S_design%*%gamma)

## set the variance of the outcome: 
outcome_sigma <- 1

### generate the outcome data 
sim_bm_outcome <- sapply(seq_along(1:I), function(i){rnorm(1, mu_outcome[i],outcome_sigma)})




## change the location if necessary: 
compiled_model <- stan_model("/home/irena/gsra/swan_sims/simulations/0511_1000_ids/joint_model.stan")



# N <- length(time_fmp)
# id <- ids
# hormones <- sim_hormones
# simulate <- 0
# 
# 
# stan_rdump(c("Q",
#                 "P",
#                 "I",
#                 "N",
#                 "id",
#                 "f_time_fmp",
#                 "hormones",
#                 "simulate"), "/home/logan/Code/irena/r/model_data.dat")

## set number of chains
no_chains = 2
## sample from the model: 
sim_out  <- sampling(compiled_model,
                     # include = TRUE,
                     sample_file=paste0('/home/irena/gsra/swan_sims/simulations/0511_1000_ids/results/rep_', taskid,'_model_samples.csv'),
                         iter = 2000,
                         warmup=1000, #BURN IN 
                         chains = no_chains,
                         seed = seed,
                         control = list(max_treedepth = 30,
                                        adapt_delta=0.95),
                         data = list(Q = Q,
                                     P = P,
                                     I = I,
                                     N = length(time_fmp),
                                     id=ids,
                                     f_time_fmp = f_time_fmp,
                                     hormones =sim_hormones,
                                     bm_outcome = sim_bm_outcome,
                                     simulate = 0
                         ))


#saveRDS(sim_out, "gsra/swan_sims/results/sim_out_500_ids.rds")
