//#######################################################
//### State-space univariate model ######################
//#######################################################

///////////////////////// DATA //////////////////////////
data {
	int<lower = 0> N;       // number of data
	real Y[N];  	// response vector
	real sigma2phi;
	real a_sigma2;
	real b_sigma2;
	real sigma2m0;
}

//////////////////// PARAMETERS /////////////////////////
parameters {
	real phi;
	real sigma2;
	real m0;
}


////////////////// MODEL ////////////////////////////////
model {

	// Likelihood     
	for (t in 2:N)
	{
		Y[t] ~ normal(m0 + phi * Y[t - 1], pow(sigma2, 0.5));
	} 
	m0 ~ normal(0, sigma2m0);
	phi ~ normal(0, sigma2phi);
	sigma2 ~ inv_gamma(a_sigma2, b_sigma2);
}

