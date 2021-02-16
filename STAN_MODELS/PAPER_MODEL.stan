// The input data is a vector 'y' of length 'N'.

data {
 int<lower=0> N; // number of states considered
 int mu_Mu;
 int<lower=0> sigma2_Mu;
 

  int<lower=1> Nwhite[N];
  int<lower=1> Nblack[N];
  int<lower=1> Nhispanic[N];

  int<lower=0> UnarmedWhite[N];
  int<lower=0> UnarmedBlack[N];
  int<lower=0> UnarmedHispanic[N];

  int<lower=0> ArmedWhite[N];
  int<lower=0> ArmedBlack[N];
  int<lower=0> ArmedHispanic[N];
}

parameters {
  vector[6] Mu;
   vector<lower=0>[6] Sigma;
   corr_matrix[6] Rho;
 
   vector[6] Theta[N];
}

model {
  /// Priors
  //Mu ~ normal(-14,4);
   //Sigma ~ cauchy(0,5);
  Mu ~ normal(mu_Mu,sigma2_Mu);
  Sigma ~ cauchy(0,5);

for(i in 1:N){
   Theta[i] ~ multi_normal_cholesky(Mu, (diag_matrix(Sigma) * cholesky_decompose(Rho)) );
}

  /// Data Modeling
for(i in 1:N){
   ArmedBlack[i]~binomial(Nblack[i],inv_logit(Theta[i,1]));
   ArmedWhite[i]~binomial(Nwhite[i],inv_logit(Theta[i,3]));
   ArmedHispanic[i]~binomial(Nhispanic[i],inv_logit(Theta[i,5]));  

   UnarmedBlack[i]~binomial(Nblack[i],inv_logit(Theta[i,2])); 
   UnarmedWhite[i]~binomial(Nwhite[i],inv_logit(Theta[i,4]));
   UnarmedHispanic[i]~binomial(Nhispanic[i],inv_logit(Theta[i,6]));
   }
}
                         
generated quantities{  /// What is this part about??? Maybe the Y~N(mu_theta,sigma_theta)...
/// Mean Quanitities
real Mu_Black_Armed; 
real Mu_White_Armed; 
real Mu_Hispanic_Armed;

real Mu_Black_Unarmed; 
real Mu_White_Unarmed; 
real Mu_Hispanic_Unarmed; 
   
real Mu_RR_Black_Armed_Versus_Unarmed; 
real Mu_RR_White_Armed_Versus_Unarmed; 
real Mu_RR_Hispanic_Armed_Versus_Unarmed;

real Mu_RR_Black_Armed_Versus_White_Armed; 
real Mu_RR_Hispanic_Armed_Versus_White_Armed; 
real Mu_RR_Hispanic_Armed_Versus_Black_Armed; 

real Mu_RR_Black_Unarmed_Versus_White_Unarmed; 
real Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed; 
real Mu_RR_Hispanic_Unarmed_Versus_Black_Unarmed; 

real Mu_RR_Black_Unarmed_Versus_White_Armed; 
real Mu_RR_Hispanic_Unarmed_Versus_White_Armed; 

/// Quanitities By County (in our case by state!!!)
vector[N]  Black_Armed; 
vector[N]  White_Armed; 
vector[N]  Hispanic_Armed;

vector[N]  Black_Unarmed; 
vector[N]  White_Unarmed; 
vector[N]  Hispanic_Unarmed; 
   
vector[N]  RR_Black_Armed_Versus_Unarmed; 
vector[N]  RR_White_Armed_Versus_Unarmed; 
vector[N]  RR_Hispanic_Armed_Versus_Unarmed;

vector[N]  RR_Black_Armed_Versus_White_Armed; 
vector[N]  RR_Hispanic_Armed_Versus_White_Armed; 
vector[N]  RR_Hispanic_Armed_Versus_Black_Armed; 

vector[N]  RR_Black_Unarmed_Versus_White_Unarmed; 
vector[N]  RR_Hispanic_Unarmed_Versus_White_Unarmed; 
vector[N]  RR_Hispanic_Unarmed_Versus_Black_Unarmed; 

vector[N]  RR_Black_Unarmed_Versus_White_Armed; 
vector[N]  RR_Hispanic_Unarmed_Versus_White_Armed; 
vector[N] log_lik_ab;
vector[N] log_lik_aw;
vector[N] log_lik_ah;
vector[N] log_lik_ub;
vector[N] log_lik_uw;
vector[N] log_lik_uh;

for (j in 1:N){
  log_lik_ab[j] = binomial_lpmf(ArmedBlack[j] |  Nblack[j], inv_logit(Theta[j,1]));
  log_lik_aw[j] = binomial_lpmf(ArmedWhite[j] |  Nwhite[j], inv_logit(Theta[j,2]));
  log_lik_ah[j] = binomial_lpmf(ArmedHispanic[j] |  Nhispanic[j], inv_logit(Theta[j,5]));
  log_lik_ub[j] = binomial_lpmf(UnarmedBlack[j] |  Nblack[j], inv_logit(Theta[j,2]));
  log_lik_uw[j] = binomial_lpmf(UnarmedWhite[j] |  Nblack[j], inv_logit(Theta[j,4]));
  log_lik_uh[j] = binomial_lpmf(UnarmedHispanic[j] |  Nhispanic[j], inv_logit(Theta[j,6]));
}

/// Calc Means
  Mu_Black_Armed = inv_logit(Mu[1]); 
  Mu_White_Armed = inv_logit(Mu[3]); 
  Mu_Hispanic_Armed = inv_logit(Mu[5]);

  Mu_Black_Unarmed = inv_logit(Mu[2]); 
  Mu_White_Unarmed = inv_logit(Mu[4]); 
  Mu_Hispanic_Unarmed = inv_logit(Mu[6]);


  Mu_RR_Black_Armed_Versus_Unarmed            =  inv_logit(Mu[1])/inv_logit(Mu[2]);
  Mu_RR_White_Armed_Versus_Unarmed            =  inv_logit(Mu[3])/inv_logit(Mu[4]);
  Mu_RR_Hispanic_Armed_Versus_Unarmed         =  inv_logit(Mu[5])/inv_logit(Mu[6]);

  Mu_RR_Black_Armed_Versus_White_Armed        =  inv_logit(Mu[1])/inv_logit(Mu[3]);
  Mu_RR_Hispanic_Armed_Versus_White_Armed     =  inv_logit(Mu[5])/inv_logit(Mu[3]);
  Mu_RR_Hispanic_Armed_Versus_Black_Armed     =  inv_logit(Mu[5])/inv_logit(Mu[1]);

  Mu_RR_Black_Unarmed_Versus_White_Unarmed    =  inv_logit(Mu[2])/inv_logit(Mu[4]);
  Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed =  inv_logit(Mu[6])/inv_logit(Mu[4]);
  Mu_RR_Hispanic_Unarmed_Versus_Black_Unarmed =  inv_logit(Mu[6])/inv_logit(Mu[2]);

  Mu_RR_Black_Unarmed_Versus_White_Armed      =  inv_logit(Mu[2])/inv_logit(Mu[3]); 
  Mu_RR_Hispanic_Unarmed_Versus_White_Armed   =  inv_logit(Mu[6])/inv_logit(Mu[3]);

  /// Calc Full Vectors
  for(i in 1:N){
  Black_Armed[i]   = inv_logit(Theta[i,1]); 
  White_Armed[i]   = inv_logit(Theta[i,3]); 
  Hispanic_Armed[i] = inv_logit(Theta[i,5]);

  Black_Unarmed[i]   = inv_logit(Theta[i,2]); 
  White_Unarmed[i]   = inv_logit(Theta[i,4]); 
  Hispanic_Unarmed[i]= inv_logit(Theta[i,6]);


  RR_Black_Armed_Versus_Unarmed[i]    =   inv_logit(Theta[i,1])/inv_logit(Theta[i,2]);
  RR_White_Armed_Versus_Unarmed[i]    =   inv_logit(Theta[i,3])/inv_logit(Theta[i,4]);
  RR_Hispanic_Armed_Versus_Unarmed[i] =   inv_logit(Theta[i,5])/inv_logit(Theta[i,6]);

  RR_Black_Armed_Versus_White_Armed[i]    =   inv_logit(Theta[i,1])/inv_logit(Theta[i,3]);
  RR_Hispanic_Armed_Versus_White_Armed[i] =   inv_logit(Theta[i,5])/inv_logit(Theta[i,3]);
  RR_Hispanic_Armed_Versus_Black_Armed[i] =   inv_logit(Theta[i,5])/inv_logit(Theta[i,1]);

  RR_Black_Unarmed_Versus_White_Unarmed[i]   =   inv_logit(Theta[i,2])/inv_logit(Theta[i,4]);
  RR_Hispanic_Unarmed_Versus_White_Unarmed[i]=   inv_logit(Theta[i,6])/inv_logit(Theta[i,4]);
  RR_Hispanic_Unarmed_Versus_Black_Unarmed[i]=   inv_logit(Theta[i,6])/inv_logit(Theta[i,2]);

  RR_Black_Unarmed_Versus_White_Armed[i]     =   inv_logit(Theta[i,2])/inv_logit(Theta[i,3]); 
  RR_Hispanic_Unarmed_Versus_White_Armed[i]  =   inv_logit(Theta[i,6])/inv_logit(Theta[i,3]);
  }
}
