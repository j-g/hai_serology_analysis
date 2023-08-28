################################################################################
# Functions to process stan model output
################################################################################

################################################################################
# Process spline output
################################################################################
process_splines<-function(out)
{
  fit<-out$fit
  N<-length(out$knots)
  
  df<-NULL
  for(j in 1:N)
  {
    y<-extract(fit)
    df_spline<-y[[paste0("a", j)]]
    
    # Set up basis and x axis
    knots<-out$knots[[j]]
    deg<-out$deg
    x<-seq(min(knots),max(knots),0.1)
    B <- t(bSpline(x, knots = knots[-c(1, length(knots))], Boundary.knots = c(knots[1], knots[length(knots)]), degree = deg, intercept = TRUE))
    
    # Sample coefficients
    spline_samp<-NULL
    for(i in 1:1000)
    {
      ind<-sample(1:nrow(df_spline), 1)
      a<-matrix(df_spline[ind,], nrow=1)
      y<-as.vector(a%*%B)
      spline_samp<-rbind(spline_samp, data.frame(x, y, i))
    }
    
    # Get CI's
    spline_mn<-spline_samp%>%group_by(x)%>%summarise(med=median(y),
                                                     lolo=quantile(y, 0.025),
                                                     lo=quantile(y, 0.25),
                                                     hi=quantile(y,0.75),
                                                     hihi=quantile(y,0.975))
    spline_mn$spline<-paste0("spline ", j)
    
    df<-rbind(df, spline_mn)
  }
  
  return(df)
}

################################################################################
# Process spline and strain output
################################################################################
process_spline_out<-function(out)
{
  #####################
  # Process spline term
  #####################
  
  spline_tsfe<-process_splines(out)
  
  ####################
  # Process Intercepts
  ####################
  
  # Get dataframe of intercepts
  dat<-out$dat
  cohort<-as.character(dat$Study[1])
  y<-extract(out$fit)
  df<-as.data.frame(y$alpha)
  df$iter<-1:nrow(df)
  df<-melt(df, id="iter")
  df$variable<-as.numeric(gsub("V","",df$variable))
  
  # Merge with actual names
  strains<-data.frame(variable=1:nlevels(dat$strain_name), strain=levels(dat$strain_name))
  df<-merge(df, strains)
  
  df<-df%>%group_by(strain)%>%summarise(med=median(value),
                                        lolo=quantile(value, 0.025),
                                        lo=quantile(value, 0.25),
                                        hi=quantile(value, 0.75),
                                        hihi=quantile(value, 0.975))
  
  
  df$year<-as.numeric( gsub(".*\\_","",df$strain) )
  
  spline_tsfe$cohort<-cohort
  df$cohort<-cohort
  
  return(list("spline_tsfe"=spline_tsfe,
              "strain"=df))
}

################################################################################
# Process parametric model
################################################################################
process_param_out<-function(out)
{
  cohort<-as.character(dat$Study[1])
  dat<-out$dat
  fit<-out$fit
  
  traces<-extract(fit)
  
  as_time<-traces$as_time
  N<-length(as_time)
  
  ##############################################################################
  # TSFE curve
  ##############################################################################
  
  as_time<-traces$as_time
  as_mag<-traces$as_mag
  
  t_max<-max(dat$TSFE)
  
  TSFE<-seq(0,t_max,0.1)
  
  df_TSFE<-NULL
  for(i in 1:1000)
  {
    ind<-sample(1:N, 1)
    
    as_time_i<-as_time[ind]
    as_mag_i<-as_mag[ind]
    
    ant_sen_c<-sqrt(2*pi*as_time_i*as_mag_i/t_max)*(0.5-pnorm(t_max, 0, as_time_i))
    y<-as_mag_i*exp(-0.5*( TSFE/as_time_i )^2)+ant_sen_c
    
    df_TSFE<-rbind(df_TSFE, data.frame(TSFE, y))
  }
  
  # Get CI's
  df_TSFE<-df_TSFE%>%group_by(TSFE)%>%summarise(med=median(y),
                                                lolo=quantile(y, 0.025),
                                                lo=quantile(y, 0.25),
                                                hi=quantile(y,0.75),
                                                hihi=quantile(y,0.975))
  
  ##############################################################################
  # TSC curve
  ##############################################################################
  
  strain0<-traces$strain0
  peak_titer<-traces$peak_titer
  t_buildup<-traces$t_buildup
  t_wane<-traces$t_wane
  
  TSC<-seq(0,max(dat$TSC),0.1)
  
  df_TSC<-NULL
  for(i in 1:1000)
  {
    ind<-sample(1:N, 1)
    
    strain0_i<-strain0[ind]
    peak_titer_i<-peak_titer[ind]
    t_wane_i<-t_wane[ind]
    t_buildup_i<-t_buildup[ind]
    
    y<-strain0_i+peak_titer_i*(1-exp(-TSC/t_buildup_i))-(TSC/t_wane_i)
    
    df_TSC<-rbind(df_TSC, data.frame(TSC, y))
  }
  
  # Get CI's
  df_TSC<-df_TSC%>%group_by(TSC)%>%summarise(med=median(y),
                                             lolo=quantile(y, 0.025),
                                             lo=quantile(y, 0.25),
                                             hi=quantile(y,0.75),
                                             hihi=quantile(y,0.975))
  
  ##############################################################################
  # AC curve
  ##############################################################################
  
  AC_a<-traces$AC_a
  AC_b<-traces$AC_b
  
  AgeCirc<-seq(0,max(dat$AgeCirc),0.1)
  
  df_AC<-NULL
  for(i in 1:1000)
  {
    ind<-sample(1:N, 1)
    
    AC_a_i<-AC_a[ind]
    AC_b_i<-AC_b[ind]
    
    y<-(AC_a_i)*(( (AgeCirc/max(dat$AgeCirc)) -AC_b_i)^2)
    
    df_AC<-rbind(df_AC, data.frame(AgeCirc, y))
  }
  
  # Get CI's
  df_AC<-df_AC%>%group_by(AgeCirc)%>%summarise(med=median(y),
                                               lolo=quantile(y, 0.025),
                                               lo=quantile(y, 0.25),
                                               hi=quantile(y,0.75),
                                               hihi=quantile(y,0.975))
  
  ##############################################################################
  # Strains
  ##############################################################################
  
  df<-as.data.frame(traces$alpha)
  df$iter<-1:nrow(df)
  df<-melt(df, id="iter")
  df$variable<-as.numeric(gsub("V","",df$variable))
  
  # Merge with actual names
  strains<-data.frame(variable=1:nlevels(dat$strain_name), strain=levels(dat$strain_name))
  df<-merge(df, strains)
  
  df_strain<-df%>%group_by(strain)%>%summarise(med=median(value),
                                               lolo=quantile(value, 0.025),
                                               lo=quantile(value, 0.25),
                                               hi=quantile(value, 0.75),
                                               hihi=quantile(value, 0.975))
  
  df_strain$year<-as.numeric( gsub(".*\\_","",df_strain$strain) )
  
  ##############################################################################
  # Raw params
  ##############################################################################
  
  y<-c(as_time, as_mag, strain0, peak_titer, t_buildup, t_wane, AC_a, AC_b)
  name<-c("as_time", "as_mag", "strain0", "peak_titer", "t_buildup", "t_wane", "AC_a", "AC_b")
  param<-rep(name, each=N)
  
  df_param<-data.frame(param, y)
  
  df_param<-df_param%>%group_by(param)%>%summarise(med=median(y),
                                                   lolo=quantile(y, 0.025),
                                                   lo=quantile(y, 0.25),
                                                   hi=quantile(y,0.75),
                                                   hihi=quantile(y,0.975))
  
  df_AC$cohort<-cohort
  df_TSFE$cohort<-cohort
  df_TSC$cohort<-cohort
  df_param$cohort<-cohort
  df_strain$cohort<-cohort
  
  return(list("df_AC"=df_AC,
              "df_TSFE"=df_TSFE,
              "df_TSC"=df_TSC,
              "df_param"=df_param,
              "df_strain"=df_strain))
}


################################################################################
# Get loo
################################################################################
get_loo<-function(out, n_cores)
{
  fit<-out$fit
  LLarray <- loo::extract_log_lik(stanfit = fit, parameter_name = "log_lik", merge_chains = FALSE)
  r_eff <- loo::relative_eff(x = exp(LLarray))
  out_loo<-loo::loo.array(LLarray, r_eff = r_eff, cores = n_cores)
  return(out_loo)
}
