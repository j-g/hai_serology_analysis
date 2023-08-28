functions {
  // Probability of measuring titer 'meas' when normally distributed
  // about mean value 'mn_pred' with sd 'sigma'. Titer is measured in
  // discrete log units, so probability is integral of cdf between meas
  // and meas+1 if in experimental range. If above upper dilution limit
  // (titer_hi) then integral is from titer_hi to +inf and if below lower
  // limit (titer_lo) then integral is from -inf to titer_lo.
  
  // Below lower limit
  real prob_lo(real titer_lo, real mn_pred, real sigma){
    return(log(normal_cdf(titer_lo+1, mn_pred, sigma)));
  }

  // Above upper limit
  real prob_hi(real titer_hi, real mn_pred, real sigma){
    real diff=titer_hi-mn_pred;
    return(log(normal_cdf(mn_pred-diff, mn_pred, sigma)));
  }
  
  // Within range, need to ensure that integral is calculated to the left
  // side of the mean to avoid rounding errors of small numbers
  real prob_in_range(real meas, real mn_pred, real sigma){
	real diff=meas-mn_pred;
    if(diff<0){
      return(log(normal_cdf(diff+1, 0, sigma)-normal_cdf(diff, 0, sigma)));
    }
    else{
      return(log(normal_cdf(-diff, 0, sigma)-normal_cdf(-(diff+1), 0, sigma)));
    }
  }
}

data {
  int <lower = 0> n_data;							// no. data points
  int <lower = 0> n_strain;							// no. strains
  int <lower = 0> n_basis1;							// no. spline basis functions for variable 1
  vector[n_data] LogTiter;							// Outcome variable, log titer
  int <lower = 0, upper = n_strain> Strain[n_data]; // Strain name (numeric)
  matrix[n_basis1, n_data] B1;						// Values of var 1 basis functions for each function at each data point
  vector[n_basis1] norms1;							// Areas under curves of var 1 basis functions
  // Likelihood calculation calls 3 different functions depending on
  // value of LogTiter. To avoid selecting the appropriate function
  // each time, data is sorted in ascending order and the indexes
  // when LogTiter exceeds the lower limit and upper limit is given.
  int <lower = 0> index_first_above_lower;								
  int <lower = 0> index_last_below_upper;
}

parameters {
  vector[n_strain] alpha;			// Strain intercepts
  real<lower = 0> sigma_y;			// Variation of outcome (measurement + individual)
  row_vector[n_basis1] a_raw1;		// Coefficients of var 1 basis functions
}

transformed parameters {
  // Normalize spline coefficients to ensure identifiability
  // by making sure integral over function is zero.
  // Smooth out curve by adding some correlation between splines.

  row_vector[n_basis1] a1;
  real scale_fac;
  
  //// Var 1
  // Smooth prior
  a1[1]=a_raw1[1];
  for (i in 2:n_basis1)
    a1[i] = a1[i-1] + a_raw1[i];
  
  // Subtract area under curve
  scale_fac=-dot_product(a1, norms1)/sum(norms1);
  a1=a1+scale_fac;
}

model {
  // Predicted log titers, y, are the sum of strain and spline terms
  vector[n_data] y;   
  y = alpha[Strain] + to_vector(a1*B1);

  // Priors
  alpha ~ normal(2,4);
  a_raw1 ~ normal(0, 1);
  sigma_y ~ exponential(1);
  
  // Compute likelihood function
  {
    vector[n_data] summands;
	
    for (i in 1:(index_first_above_lower-1)){
      summands[i] = prob_lo(LogTiter[i], y[i], sigma_y);
    }
	
	for (i in index_first_above_lower:index_last_below_upper){
      summands[i] = prob_in_range(LogTiter[i], y[i], sigma_y);
    }
	
	if(index_last_below_upper+1<n_data){
	  for (i in (index_last_below_upper+1):n_data){
        summands[i] = prob_hi(LogTiter[i], y[i], sigma_y);
      }
	}
    target += sum(summands);
  }
}

generated quantities {
  // Quantities used for leave-one-out-cross-validation
  vector[n_data] log_lik;
  
  {
    vector[n_data] y;
  
    y = alpha[Strain] + to_vector(a1*B1);

	for (i in 1:(index_first_above_lower-1)){
      log_lik[i] = prob_lo(LogTiter[i], y[i], sigma_y);
    }
	
	for (i in index_first_above_lower:index_last_below_upper){
      log_lik[i] = prob_in_range(LogTiter[i], y[i], sigma_y);
    }
	
	for (i in (index_last_below_upper+1):n_data){
      log_lik[i] = prob_hi(LogTiter[i], y[i], sigma_y);
    }
	
  }

}
