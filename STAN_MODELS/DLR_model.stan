//#######################################################
//### State-space univariate model ######################
//#######################################################

///////////////////////// DATA //////////////////////////
data {
	int<lower = 0> N;       // number of data
	real Y[N];  	// response vector
	real X[N];  	// response vector
	real sigma_coef[3];
	real sigma2Y_par[2];
}

//////////////////// PARAMETERS /////////////////////////
parameters {
	real sigma2Y;
	real m0Y;
	real beta;
	real gamma;
}


////////////////// MODEL ////////////////////////////////
model {

	// Likelihood     
	for (t in 2:N)
	{
		Y[t] ~ normal(m0Y + beta * X[t] + gamma * Y[t - 1], pow(sigma2Y, 0.5));
	} 
	
	m0Y ~ normal(0, sigma_coef[1]);
	beta ~ normal(0, sigma_coef[2]);
	gamma ~ normal(0, sigma_coef[3]);
	sigma2Y ~ inv_gamma(sigma2Y_par[1], sigma2Y_par[2]);
}

