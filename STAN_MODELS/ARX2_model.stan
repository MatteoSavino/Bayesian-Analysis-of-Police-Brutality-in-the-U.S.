//
//#######################################################
//### State-space univariate model ######################
//#######################################################

///////////////////////// DATA //////////////////////////
data {
	int<lower = 0> N;       // number of data
	real Y[N];  	// response vector
	real X[N];  	// response vector
	real Z[N];  	// response vector
	real sigma_coef[5];
	real sigma2Y_par[2];
}

//////////////////// PARAMETERS /////////////////////////
parameters {
	real sigma2Y;
	real m0Y;
	real beta1;
	real beta2;
	real gamma1;
	real gamma2;
}


////////////////// MODEL ////////////////////////////////
model {

	// Likelihood     
	for (t in 3:N)
	{
		Y[t] ~ normal(m0Y + beta1 * X[t] + beta2 * Z[t] + gamma1 * Y[t - 1] + gamma2 * Y[t - 2], pow(sigma2Y, 0.5));
	} 
	
	m0Y ~ normal(0, sigma_coef[1]);
	beta1 ~ normal(0, sigma_coef[2]);
	beta2 ~ normal(0, sigma_coef[3]);
	gamma1 ~ normal(0, sigma_coef[4]);
	gamma2 ~ normal(0, sigma_coef[5]);
	sigma2Y ~ inv_gamma(sigma2Y_par[1], sigma2Y_par[2]);
}
