// this is to get the B_is and S_is into a shape that can be used in the second submodel
functions {
  matrix to_triangular(vector x, int K) {
    // could check rows(y) = K * (K + 1) / 2
    matrix[K, K] y;
    int pos = 1;
    for (i in 1: K) {
      for (j in 1:i) {
        y[i, j] = x[pos];
        pos += 1;
      }
      for (j in (i+1):K) {
        y[i, j] = 0;
      }
    }
    return y;
  }
}

data {
  int<lower=2> P; // df for basis expansion on time_fmp
  int<lower=1> N; // no of data points 
  int<lower=1> Q; //  no of hormones
  int<lower=1> I; //  no of subjects
  int<lower=1> id[N]; //array of subject ids (length N in order to do the longitudinal estimation)
  // data for the longitudinal submodel:
  matrix[N,P] f_time_fmp; // basis expansion on time 
  matrix[N,Q] hormones; // this should be X_i at time t (Q-length vector of hormone values)
  
  //data for the outcome submodel:
  vector[I] bm_outcome;
  // Whether or not to evaluate the likelihood
  int<lower = 0, upper = 1> simulate; 
}

parameters {
  //for the longitudinal submodel 
  vector[P] beta[Q]; // mean on the B_ijs 
  // vector[P] beta_raw[Q];
  matrix[P,Q] B[I]; // Q vectors of growth curves
  // parameters for Sigma 
  vector<lower=0>[Q] S_tau; // prior scale (for covariance matrix for the Bi's)
  corr_matrix[Q] S_Omega;
  cholesky_factor_corr[P] L_Omega[Q];  // prior correlation (for the covariance matrix of Bijs)
  vector<lower=0, upper=pi()/2>[P] tau_unif[Q];
  //parameters for Si 
  vector[Q] log_Q_Si[I]; // prior scale for Si
  real<lower=0, upper=1> t_r_beta[I];
  real<lower=0> hyper_alpha; // hyperparameter on the Beta prior
  real<lower=0> hyper_beta; //hyperparameter on the Beta prior
  vector[Q] hyper_mu; //hyperparameter on the log-Normal
  vector<lower=0, upper=pi()/2>[Q] tau_hyper_sigma; //hyperparamter on the log-Normal
  
 //for the outcome submodel 
   real a0; //intercept for the outcome 
  matrix[P,Q] alpha; //coefficients for the mean 
  vector[Q + choose(Q, 2)] gamma_flat; //coefficients for the (co)variances, size of lower triangle of S
  real<lower=0> outcome_sigma; //parameter for the variance of the outcomes

}

transformed parameters{
   cov_matrix[Q] S[I]; //subject specific vcov 
   vector<lower=0>[Q] Q_Si[I]; // standard deviations 
   matrix[Q,Q] R_Si[I]; // corr matrix for Si
   real<lower=-1, upper=1> r_beta[I]; // 
   cov_matrix[P] Sigma[Q]; // vcov on the B_ijs  
   vector<lower=0>[P] tau[Q]; // prior scale (for covariance matrix for the Bi's)
   vector<lower=0>[Q] hyper_sigma; //hyperparamter on the log-Normal
   matrix[Q,Q] gamma;
   
   // get the gamma parameters as a matrix 
    gamma = to_triangular(gamma_flat, Q);
    
  // get back the population variance-covariance of the hormone trajectories 
  for(q in 1:Q){
     tau[q] = 2.5 * tan(tau_unif[q]);
     hyper_sigma[q] =tan(tau_hyper_sigma[q]);
     Sigma[q] =  diag_pre_multiply(tau[q], L_Omega[q]) * diag_pre_multiply(tau[q], L_Omega[q])';
   }  
   
  for(i in 1:I){
    r_beta[i] = 2*(t_r_beta[i])-1; //we want this to range between [-1,1] 
    //eventually will code this better into an upper triangular matrix for Q>2 but I think this works for just a 2x2 matrix 
    R_Si[i][2,1] =r_beta[i]; //offdiagonals of the correlation matrix should be equal 
    R_Si[i][1,2] = r_beta[i]; //offdiagonals of the correlation matrix should be equal 
    R_Si[i][1,1] = 1; //diagonals on the correlation matrix are 1 
    R_Si[i][2,2] = 1; //diagonals on the correlation matrix are 1 
    Q_Si[i] = exp(log_Q_Si[i]); // get the qi values from the log-normal
    S[i] = quad_form_diag(R_Si[i], Q_Si[i]); //get back the S_i's after decomposing
  }
}

model {
  
  //hyperparameters on the S_is 
  hyper_alpha ~ exponential(0.1);
  hyper_beta ~ exponential(0.1);
  hyper_mu ~ normal(0, 10);
  
   //prior on variance of the outcome model 
  outcome_sigma ~ cauchy(0,2.5); //weak prior on the sigma for outcome model
  
    //outcome submodel
  to_vector(alpha) ~ normal(0, 10);
  to_vector(gamma_flat) ~ normal(0, 10); 
  //set up the variance parameters for each beta: 
  for(q in 1:Q){
      // rather than put an inverse wishart on Sigma, follow stan recommendations for LKJ: 
      L_Omega[q] ~ lkj_corr_cholesky(1);
      // draw betas from a nmultinormal distribution
      beta[q] ~ normal(0,10); //put a WEAKLY informative prior
      //   beta_raw[p,q] ~ std_normal();
  }
  
  for(i in 1:I) {
     // draw values for the individual Si's: 
    t_r_beta[i] ~ beta(hyper_alpha,hyper_beta); // prior on the correlation matrix
   
   // scale from a lognormal distribution
    for(q in 1:Q){
       log_Q_Si[i][q] ~ normal(hyper_mu[q], hyper_sigma[q]); 
      //draw individual mean parameters from population distribution
      B[i][,q] ~ multi_normal(beta[q], Sigma[q]); //instaed of sigma[q]
    }
  }
  
    
  S_Omega ~ lkj_corr(1);
  S_tau ~ cauchy(0,2.5);
  // this estimates the parameters using the time fmp and hormone data 
  {
  if(simulate == 0){ 
       vector[Q] mu[N]; // for each person, vector of two hormone means
       real bm_outcome_mu[I];
      for(n in 1:N){
        for(q in 1:Q){ // for each hormone
           mu[n][q] = dot_product(B[id[n]][,q], f_time_fmp[n]);  //to get the means of the bi's
        }
         hormones[n] ~ multi_normal(mu[n], S[id[n]]);
      }
     for(i in 1:I){
      bm_outcome_mu[i] = a0 + sum(B[i] .* alpha) + sum(S[i] .* gamma);
      bm_outcome[i] ~ normal(bm_outcome_mu[i], outcome_sigma);
     }
   }
  } 
}



// generated quantities {
//   // uncomment if you want to simulate data! 
//   matrix[N,Q] sim_hormones;
//   vector[Q] mu[N]; 
//   vector[I] sim_outcome;
//   vector[I] outcome_mu;
//    // uncomment if you want to simulate data! 
//    //simulate hormone data instead of using the data to generate it 
//   if(simulate == 1){
//     for(n in 1:N){
//        for(q in 1:Q){ // for each hormone
//            mu[n][q] = dot_product(B[id[n]][,q], f_time_fmp[n]);  //to get the means of the bi's
//       }
//       sim_hormones[n] = to_row_vector(multi_normal_rng(mu[n], S[id[n]]));
//     }
//     for(i in 1:I){
//       outcome_mu[i] = sum(B[i] .* alpha) + sum(S[i] .* gamma);
//       sim_outcome[i] = normal_rng(outcome_mu[i], outcome_sigma);
//     }
//   }
// 
// }






