/// POISSON MODEL ///

///////////////////////// DATA /////////////////////////////////////
data {
  int<lower=0> N; // number of states considered
  int<lower=0> p_fix; // number of covariate (fixed effects)
  int<lower = 0> p_ran;    // number of covariate (random effects)
  int<lower = 0> group[N];    

	matrix[N,p_fix] X;   	// design matrix (fixed effects)
	matrix[N,p_ran] W;   	// design matrix (random effects)
	
	real mu_beta;
  real a_sigma2_beta;
	real b_sigma2_beta;

	real mu_gamma;
	real a_sigma2_gam;
	real b_sigma2_gam;

  //int<lower=1> N_race[N];
  int<lower=0> Y[N];
  
  // offset
  real offset[N];


}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
  
   vector[p_fix] beta;
   vector[p_fix] sigma2_beta;

   vector[p_ran] gamma;
   vector[p_ran] sigma2_gamma;
   
   //vector[p_ran] u;

}

//////////////////// TRANSFORMED PARAMETERS ///////////////////////
transformed parameters{
  
 vector[N] mu; // N x 1
 for(i in 1:N){
   
    mu[i] = exp(row(W,i) * gamma + row(X,i) * beta + offset[i]); 
    
    
	}
	
}

////////////////// MODEL ////////////////////////
model {

  /// LIKELIHOOD
for(i in 1:N){
		Y[i] ~ poisson(mu[i]);  
   }
   
  /// PRIORS
	for (i in 1:p_fix) {
	 	beta[i] ~ normal(mu_beta, pow(sigma2_beta[i], 0.5));
		sigma2_beta[i] ~ inv_gamma(a_sigma2_beta, b_sigma2_beta);
	}
	
   for(i in 1:p_ran){
     gamma[i] ~ normal(mu_gamma,pow(sigma2_gamma[i], 0.5));
     sigma2_gamma[i] ~ inv_gamma(a_sigma2_gam,b_sigma2_gam);
  }
  
  //for(i in 1:p_ran){
  //  u[i] ~ normal(0,1);
  //}
  
}

///////////////////////// GENERATED QUANTITIES ////////////////////
generated quantities{
  
vector[N] log_lik;

for (j in 1:N){
  log_lik[j] = poisson_lpmf(Y[j] | mu[j]);
  }

}

