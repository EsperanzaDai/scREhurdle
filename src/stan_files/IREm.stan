data {
  int<lower=0> G;               // number of genes
  int<lower=0> N;               // number of cells
  int<lower=1> M;               // number of columns in model matrix X + 1 to include random effects
  int<lower=0> y[G,N];          // expression matrix
  matrix[N, M-1] X;             // model matrix
  vector<lower=0>[M] prior_s;   // scale for cauchy priors (betas/zetas)
}// end data

parameters {
  matrix[M,G] beta_L;           // logistic betas: beta_L[M,] = zeta_L
  vector[N] omega;              // independent random effects        
  real<lower=0> sigma2;         // variance of random effects
}// end parameters

model {
  matrix[N,G] X_beta_L;
  matrix[N,M] X_full;
  X_full = append_col(X,omega);
  
  X_beta_L = X_full*beta_L;
  
  for(g in 1:G){
    for(i in 1:N){
      if (y[g,i] == 0) 
        1 ~ bernoulli_logit(X_beta_L[i,g]);
      else{
        0 ~ bernoulli_logit(X_beta_L[i,g]);
      }
    }// end cell loop
  }// end gene loop
    
  sigma2 ~ inv_gamma(1,1);
  omega ~ normal(0,sqrt(sigma2));

  for(j in 1:M){
    beta_L[j,] ~ cauchy(0,prior_s[j]);
  }
    
}// end model

generated quantities {
  row_vector[G] beta_L1;        // beta_L corresponding to treatment effect
  row_vector[G] zeta_L;         // direct access zeta_L (beta_L[M,] = zeta_L)
  beta_L1 = beta_L[2,];
  zeta_L = beta_L[M,];
}// end gq


