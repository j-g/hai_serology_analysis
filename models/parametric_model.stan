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
  int <lower = 0> n_data;                             // no. data points
  int <lower = 0> n_strain;                           // no. strains
  vector[n_data] LogTiter;                            // Outcome variable LogTiter
  vector[n_data] TSFE;                                // ImmAgeCirc
  vector[n_data] AgeCirc;                             // AgeCirc normalized to max value
  vector[n_data] TSC;                                 // TSC
  int <lower = 0, upper = n_strain> Strain[n_data];   // Strain name (numeric)
  real ImmAgeCirc_max;                                // Max value of ImmAgeCirc
  // Likelihood calculation calls 3 different functions depending on
  // value of LogTiter. To avoid selecting the appropriate function
  // each time, data is sorted in ascending order and the indexes
  // when LogTiter exceeds the lower limit and upper limit is given.
  int <lower = 0> index_first_above_lower;								
  int <lower = 0> index_last_below_upper;
}

parameters {
  // TSC function parameters (see manuscript methods)
  real<lower = 0> strain0;
  real<lower = 0> peak_titer;
  real<lower = 0> t_buildup;
  real<lower = 0> t_wane;
  
  // ImmAgeCirc parameters
  real<lower = 0> as_mag;
  real<lower = 0> as_time;
  
  // AgeCirc parameters
  real<lower = 0> AC_a;
  real<lower = 0> AC_b;
  
  // Variation of outcome (measurement + individual)
  real<lower = 0> sigma_y;
  
  // Strain intercepts (now assumed random with sd sigma_a)
  vector[n_strain] alpha;
  real<lower = 0> sigma_a;
}

transformed parameters { 
  vector[n_data] Y_ImmAgeCirc;
  vector[n_data] Y_TSC;
  vector[n_data] Y_AgeCirc;
  real ant_sen_c;
  
  // ImmAgeCirc curve
  ant_sen_c=(sqrt(2*pi())*as_time*as_mag/ImmAgeCirc_max)*(0.5-normal_cdf(ImmAgeCirc_max, 0, as_time));
  Y_ImmAgeCirc=as_mag*exp(-0.5*square( TSFE/as_time ))+ant_sen_c;

  // AgeCirc curve
  Y_AgeCirc=(AC_a)*(square(AgeCirc-AC_b));

  // TSC curve
  Y_TSC=strain0+peak_titer*(1-exp(-TSC/t_buildup))-(TSC/t_wane);
}

model {
  vector[n_data] y;
  y = Y_ImmAgeCirc + Y_AgeCirc + Y_TSC + alpha[Strain];

  // Priors
  sigma_y ~ exponential(1);

  sigma_a ~ exponential(1);
  alpha ~ normal(0, sigma_a);

  strain0~normal(2,2);
  peak_titer~normal(5,5);
  t_buildup~normal(5,5);
  t_wane~normal(5,5);
  
  as_mag~normal(0,1);
  as_time~normal(0,10);
  
  AC_a~normal(0,1);
  AC_b~normal(0,1);
  
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
  vector[n_data] log_lik;
  
  {
    vector[n_data] y;
  
    y = Y_ImmAgeCirc + Y_TSC + Y_AgeCirc + alpha[Strain];
  
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
