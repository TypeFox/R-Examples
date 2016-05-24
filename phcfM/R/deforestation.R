#===================================================================
# deforestation.R
#
# deforestation() is a R function sampling from the posterior
# distribution of a logistic regression model with variable
# time-interval between censuses. The function uses embedded C++ code
# in Scythe
#
#===================================================================
#
# Original code by Ghislain Vieilledent, March 2012
# CIRAD UR BSEF
# ghislain.vieilledent@cirad.fr / ghislainv@gmail.com
#
#===================================================================
# 
# This software is distributed under the terms of the GNU GENERAL
# PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
# file for more information.
#
# Copyright (C) 2011 Andrew D. Martin and Kevin M. Quinn
#
#===================================================================
#
# Revisions: 
# - G. Vieilledent, March 2012
#
#===================================================================


deforestation <- function (formula, interval=1, data, burnin=1000, mcmc=1000,
                           thin=1, verbose=1, seed=NA, tune=1,
                           beta.start=NA, mubeta=0, Vbeta=1.0E6) {
    
  #========
  # Basic checks
  #========
  check.mcmc.parameters(burnin, mcmc, thin)
  check.verbose(verbose)

  #========
  # Seed
  #========
  seed <- form.seeds(seed)
  
  #======== 
  # Form response and model matrices
  #========
  mf <- model.frame(formula=formula,data=data)
  X <- model.matrix(attr(mf,"terms"),data=mf)
  Y <- model.response(mf)
  check.Y.Binomial(Y)
  Int <- check.interval(interval,nrow(X))

  #======== 
  # Model parameters
  #========
  nobs <- nrow(X)
  np <- ncol(X)
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin

  #========
  # Form and check starting parameters
  #========
  beta.start <- form.beta.start(formula,data,beta.start,np,family="binomial",defaults=NA)
     
  #========
  # Form priors
  #========
  mvn.prior <- form.mvn.prior(mubeta,Vbeta,np)
  mubeta <- mvn.prior[[1]]
  Vbeta <- mvn.prior[[2]]

  #========
  # Parameters to save
  #========
  beta_vect <- rep(c(beta.start),each=nsamp)
  Deviance <- rep(0,nsamp)
  
  #========
  # Tuning
  #========
  tune <- check.tune(tune)
  VCV <- vcov(glm(formula=formula, data=data, family="binomial"))

  #========
  # call C++ code to draw sample
  #========
  Sample <- .C("deforestation",
               #= Constants and data
               ngibbs=as.integer(ngibbs), nthin=as.integer(nthin), nburn=as.integer(nburn),## Number of iterations, burning and samples
               nobs=as.integer(nobs), ## Constants
               np=as.integer(np), ## Constants
               Y_vect=as.integer(c(Y)), ## Response variable
               X_vect=as.double(c(X)), ## Covariates
               Int_vect=as.double(c(Int)), ## Covariates
               #= Object for proposal in metropolis
               tune_scalar=as.double(tune),
               VCV_vect=as.double(c(VCV)),
               #= Parameters to save
               beta_vect.nonconst=as.double(beta_vect), ## Fixed parameters of the regression
               #= Defining priors
               mubeta_vect=as.double(c(mubeta)), Vbeta_vect=as.double(c(Vbeta)),
               #= Diagnostic
               Deviance.nonconst=as.double(Deviance),
               #= Seeds
               seed=as.integer(seed), 
               #= Verbose
               verbose=as.integer(verbose),
               PACKAGE="phcfM")
 
  #= Matrix of MCMC samples
  Matrix <- matrix(NA,nrow=nsamp,ncol=np)
  names.fixed <- paste("beta.",colnames(X),sep="")
  colnames(Matrix) <- names.fixed

  #= Filling-in the matrix
  Matrix[,c(1:np)] <- matrix(Sample[[11]],ncol=np)

  #= Transform Sample list in an MCMC object
  MCMC <- mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)
  
  #= Output
  return (list(mcmc=MCMC,deviance=mean(Sample[[14]]),tune=Sample[[9]]))

}

#===================================================================
# END
#===================================================================



