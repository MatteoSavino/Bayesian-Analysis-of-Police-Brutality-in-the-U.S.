//#######################################################
//### ARMA(1,1) univariate model   ######################
//#######################################################

///////////////////////// DATA //////////////////////////
data {
	int<lower = 0> N;       // number of data
	real Y[N];  	// response vector
	//real sigma2phi;
	real a_sigma2;
	real b_sigma2;
	real sigma2m0;
}

//////////////////// PARAMETERS /////////////////////////
parameters {
//	real phi;
	real theta; // MA coeff!!!
	real sigma2;
	real m0;
}

transformed parameters {
  real epsilon[N];       // error terms
  epsilon[1] = Y[1] - m0;
  for (t in 2:N)
    epsilon[t] = ( Y[t] - m0
                    - theta * epsilon[t - 1]
                    //- phi * Y[t-1]
                    );
}

////////////////// MODEL ////////////////////////////////

model {

	// Likelihood     
	for (t in 2:N)
	{
		Y[t] ~ normal(m0+ theta*epsilon[t-1], pow(sigma2, 0.5));
		//Y[t] ~ normal(m0 + phi*Y[t-1] + theta*epsilon[t-1], pow(sigma2, 0.5));

	} 
	m0 ~ normal(0, sigma2m0);
	//phi ~ normal(0, sigma2phi);
	theta ~ cauchy(0,2.5);
	sigma2 ~ inv_gamma(a_sigma2, b_sigma2);
}
