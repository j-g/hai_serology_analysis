################################################################################
# Functions to fit 3 models:
# (1) One spline model
# (2) Two spline model
# (3) Parameterised model
################################################################################

################################################################################
# 1 spline
################################################################################
fit_mod_1_spline<-function(dat, n_knot, k_space,  deg, niter, cores, var)
{
  ##############################################################################
  # Prep data
  ##############################################################################
  
  # No strains from future and nothing from before birth
  dat<-dat%>%filter(AgeCirc>=0 & TSC>=0)
  dat$Strain<-factor(dat$Strain)
  
  # Strains names need to be numerical for stan, so take a note here to remember
  # the names
  dat$strain_name<-dat$Strain
  dat$Strain<-as.numeric(dat$Strain)
  
  # Organise data in ascending order and split into measurements (1) below the
  # min titer limit of 1 (2) within experimental range and (3) above the max
  # titer of 8.
  dat<-dat[order(dat$LogTiter),]
  index_first_above_lower<-min(which(dat$LogTiter>0))
  index_last_below_upper<-min(which(dat$LogTiter>7))-1
  # However, Aus data has no upper limit as titrations were performed all the
  # way
  if(grepl("Aus", aus_97$Study[1])){
    index_last_below_upper<-nrow(dat)
  }
  
  # Sizes of objects
  n_data <- nrow(dat)
  n_strain <- nlevels(as.factor(dat$Strain))
  
  ##############################################################################
  # Prep spline basis
  ##############################################################################
  
  # Either have knots in even intervals, or in intervals proportional to the
  # data distribution
  if(k_space=="even"){knots<-seq(from=min(dat[,var]), to=max(dat[,var]), length=n_knot)}
  if(k_space=="prop"){knots<-unname(quantile(dat[,var],probs=seq(from=0, to=1, length.out = n_knot)))}
  
  # Get b-spline basis functions
  B <- t(bSpline(dat[,var], knots = knots[-c(1, length(knots))], Boundary.knots = c(knots[1], knots[length(knots)]), degree = deg, intercept = TRUE))
  
  # Find integrals of splines, used to ensure identifiability
  range<-c(min(dat[,var]), max(dat[,var]))
  bsMat <- ibs(range, knots = knots[-c(1, length(knots))], Boundary.knots = c(knots[1], knots[length(knots)]), degree = deg, intercept = TRUE)
  norms<-as.vector(bsMat[nrow(bsMat),])
  
  n_basis<-nrow(B)
  
  ##############################################################################
  # Fit model
  ##############################################################################
  
  # Stan options
  rstan_options(auto_write = TRUE)
  options(mc.cores = cores)
  
  fit <- stan(file="one_spline_model.stan",
              iter=niter,
              data=list(n_data=n_data,
                        n_strain=n_strain,
                        n_basis1=n_basis,
                        LogTiter=dat$LogTiter,
                        Strain=dat$Strain,
                        B1=B,
                        norms1=norms,
                        index_first_above_lower=index_first_above_lower,
                        index_last_below_upper=index_last_below_upper))
  
  ##############################################################################
  # Output
  ##############################################################################
  
  out<-list()
  out$dat<-dat
  out$fit<-fit
  out$deg<-deg
  out$knots<-list(knots)
  out$extr<-extract(fit)
  out$n_spline<-1
  
  return(out)
}


################################################################################
# 2 splines
################################################################################
fit_mod_2_spline<-function(dat, n_knot, k_space,  deg, niter, cores, vars)
{
  ##############################################################################
  # Prep data
  ##############################################################################
  
  # No strains from future and nothing from before birth
  dat<-dat%>%filter(AgeCirc>=0 & TSC>=0)
  dat$Strain<-factor(dat$Strain)
  
  # Strains names need to be numerical for stan, so take a note here to remember
  # the names
  dat$strain_name<-dat$Strain
  dat$Strain<-as.numeric(dat$Strain)
  
  # Organise data in ascending order and split into measurements (1) below the
  # min titer limit of 1 (2) within experimental range and (3) above the max
  # titer of 8.
  dat<-dat[order(dat$LogTiter),]
  index_first_above_lower<-min(which(dat$LogTiter>0))
  index_last_below_upper<-min(which(dat$LogTiter>7))-1
  # However, Aus data has no upper limit as titrations were performed all the
  # way
  if(grepl("Aus", aus_97$Study[1])){
    index_last_below_upper<-nrow(dat)
  }
  
  # Sizes of objects
  n_data <- nrow(dat)
  n_strain <- nlevels(as.factor(dat$Strain))
  
  ##############################################################################
  # Prep spline basis
  ##############################################################################
  B_list<-list()
  norms_list<-list()
  n_basis_list<-list()
  knots_list<-list()
  
  for(i in 1:2)
  {
    var<-vars[i]
    
    # Either have knots in even intervals, or in intervals proportional to the
    # data distribution
    if(k_space=="even"){knots<-seq(from=min(dat[,var]), to=max(dat[,var]), length=n_knot)}
    if(k_space=="prop"){knots<-unname(quantile(dat[,var],probs=seq(from=0, to=1, length.out = n_knot)))}
    
    # Get b-spline basis functions
    B <- t(bSpline(dat[,var], knots = knots[-c(1, length(knots))], Boundary.knots = c(knots[1], knots[length(knots)]), degree = deg, intercept = TRUE))
    
    # Find integrals of splines, used to ensure identifiability
    range<-c(min(dat[,var]), max(dat[,var]))
    bsMat <- ibs(range, knots = knots[-c(1, length(knots))], Boundary.knots = c(knots[1], knots[length(knots)]), degree = deg, intercept = TRUE)
    norms<-as.vector(bsMat[nrow(bsMat),])
    
    B_list[[i]]<-B
    norms_list[[i]]<-norms
    n_basis_list[[i]] <- nrow(B)
    knots_list[[i]]<-knots
  }
  
  ##############################################################################
  # Fit model
  ##############################################################################
  
  # Stan options
  rstan_options(auto_write = TRUE)
  options(mc.cores = cores)
  
  fit <- stan(file="two_spline_model.stan",
              iter=niter,
              data=list(n_data=n_data,
                        n_strain=n_strain,
                        n_basis1=n_basis_list[[1]],
                        n_basis2=n_basis_list[[1]],
                        LogTiter=dat$LogTiter,
                        Strain=dat$Strain,
                        B1=B_list[[1]],
                        B2=B_list[[2]],
                        norms1=norms_list[[1]],
                        norms2=norms_list[[2]],
                        index_first_above_lower=index_first_above_lower,
                        index_last_below_upper=index_last_below_upper))
  
  ##############################################################################
  # Output
  ##############################################################################
  
  out<-list()
  out$dat<-dat
  out$fit<-fit
  out$deg<-deg
  out$knots<-knots_list
  out$extr<-extract(fit)
  out$n_spline<-2
  
  return(out)
}


################################################################################
# Parametric model
################################################################################
fit_mod_parametric<-function(dat, niter, cores)
{
  ##############################################################################
  # Prep data
  ##############################################################################
  
  # No strains from future and nothing from before birth
  dat<-dat%>%filter(AgeCirc>=0 & TSC>=0)
  dat$Strain<-factor(dat$Strain)
  
  # Strains names need to be numerical for stan, so take a note here to remember
  # the names
  dat$strain_name<-dat$Strain
  dat$Strain<-as.numeric(dat$Strain)
  
  # Organise data in ascending order and split into measurements (1) below the
  # min titer limit of 1 (2) within experimental range and (3) above the max
  # titer of 8.
  dat<-dat[order(dat$LogTiter),]
  index_first_above_lower<-min(which(dat$LogTiter>0))
  index_last_below_upper<-min(which(dat$LogTiter>7))-1
  # However, Aus data has no upper limit as titrations were performed all the
  # way
  if(grepl("Aus", aus_97$Study[1])){
    index_last_below_upper<-nrow(dat)
  }
  
  # Sizes of objects
  n_data <- nrow(dat)
  n_strain <- nlevels(as.factor(dat$Strain))
  
  ##############################################################################
  # Fit model
  ##############################################################################
  
  # Stan options
  rstan_options(auto_write = TRUE)
  options(mc.cores = cores)
  
  fit <- stan(file="parametric_model.stan",
              iter=niter,
              data=list(n_data=n_data,
                        n_strain=n_strain,
                        LogTiter=dat$LogTiter,
                        TSFE=dat$TSFE,
                        AgeCirc=dat$AgeCirc/max(dat$AgeCirc),
                        TSC=dat$TSC,
                        Strain=dat$Strain,
                        ImmAgeCirc_max=max(dat$TSFE),
                        index_first_above_lower=index_first_above_lower,
                        index_last_below_upper=index_last_below_upper))
  
  ##############################################################################
  # Output
  ##############################################################################
  
  out<-list()
  out$dat<-dat
  out$fit<-fit
  out$extr<-extract(fit)
  
  return(out)
}