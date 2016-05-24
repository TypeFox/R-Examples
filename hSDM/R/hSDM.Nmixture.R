####################################################################
##
## hSDM.Nmixture.R
##
####################################################################
##
## Original code by Ghislain Vieilledent, November 2013
## CIRAD UR B&SEF
## ghislain.vieilledent@cirad.fr / ghislainv@gmail.com
##
####################################################################
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3. See the package LICENSE file for more
## information.
##
## Copyright (C) 2011 Ghislain Vieilledent
## 
####################################################################


hSDM.Nmixture <- function (# Observations
                           counts, observability,
                           site, data.observability,
                           # Habitat
                           suitability, data.suitability,
                           # Predictions
                           suitability.pred=NULL,
                           # Chains
                           burnin=5000, mcmc=10000, thin=10,
                           # Starting values
                           beta.start,
                           gamma.start,
                           # Priors
                           mubeta=0, Vbeta=1.0E6,
                           mugamma=0, Vgamma=1.0E6,
                           # Various
                           seed=1234, verbose=1,
                           save.p=0, save.N=0)

{   
  #========
  # Basic checks
  #========
  check.mcmc.parameters(burnin, mcmc, thin)
  check.verbose(verbose)
  check.save.p(save.p)
  check.save.N(save.N)
   
  #======== 
  # Form response, covariate matrices and model parameters
  #========

  #= Response
  Y <- counts
  nobs <- length(Y)
  #= Suitability
  mf.suit <- model.frame(formula=suitability,data=data.suitability)
  X <- model.matrix(attr(mf.suit,"terms"),data=mf.suit)
  #= Observability
  mf.obs <- model.frame(formula=observability,data=data.observability)
  W <- model.matrix(attr(mf.obs,"terms"),data=mf.obs)
  #= Spatial entity
  Levels.site <- sort(unique(site))
  nsite <- length(Levels.site)
  sites <- as.numeric(as.factor(site))

  #= Predictions
  if (is.null(suitability.pred)) {
      X.pred <- X
      npred <- nsite
  }
  if (!is.null(suitability.pred)) {
      mf.pred <- model.frame(formula=suitability,data=suitability.pred)
      X.pred <- model.matrix(attr(mf.pred,"terms"),data=mf.pred)
      npred <- nrow(X.pred)
  }
  #= Model parameters
  np <- ncol(X)
  nq <- ncol(W)
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin

  #========== 
  # Check data
  #==========
  check.Y.poisson(Y)
  check.X(X,nsite) # X must be of dim (nsite x np) for the N-mixture model
  check.W(W,nobs)
  check.sites(sites,nobs)

  #========
  # Initial starting values for M-H
  #========
  beta.start <- form.beta.start(beta.start,np)
  gamma.start <- form.gamma.start(gamma.start,nq)
  # For N, we compute the MAX of the observations on each site
  N.start <- rep(0,nsite)
  Levels.sites <- sort(unique(sites))
  for (i in 1:length(Levels.sites)) {
      N.start[Levels.sites[i]] <- max(Y[sites==Levels.sites[i]]) # ! Levels.sites here
  }
  
  #========
  # Form and check priors
  #========
  mubeta <- check.mubeta(mubeta,np)
  Vbeta <- check.Vbeta(Vbeta,np)
  mugamma <- check.mugamma(mugamma,nq)
  Vgamma <- check.Vgamma(Vgamma,nq)

  #========
  # Parameters to save
  #========
  beta <- rep(beta.start,nsamp)
  gamma <- rep(gamma.start,nsamp)
  lambda_latent <- rep(0,nsite)
  delta_latent <- rep(0,nobs)
  if (save.p==0) {lambda_pred <- rep(0,npred)}
  if (save.p==1) {lambda_pred <- rep(0,npred*nsamp)}
  Deviance <- rep(0,nsamp)
  if (save.N==0) {N_pred <- rep(0,nsite)}
  if (save.N==1) {N_pred <- rep(0,nsite*nsamp)}

  #========
  # call C++ code to draw sample
  #========
  Sample <- .C("hSDM_Nmixture",
               #= Constants and data
               ngibbs=as.integer(ngibbs), nthin=as.integer(nthin), nburn=as.integer(nburn), ## Number of iterations, burning and samples
               nobs=as.integer(nobs),
               nsite=as.integer(nsite),
               np=as.integer(np),
               nq=as.integer(nq),
               Y_vect=as.integer(c(Y)),
               W_vect=as.double(c(W)),
               X_vect=as.double(c(X)),
               #= Spatial sites
               C_vect=as.integer(c(sites)-1), # Sites range is 1,...,nsite in R. Must start at 0 for C. Don't forget the "-1" term. 
               #= Predictions
               npred=as.integer(npred),
               X_pred_vect=as.double(c(X.pred)),
               #= Starting values for M-H
               beta_start=as.double(c(beta.start)),
               gamma_start=as.double(c(gamma.start)),
               N_start=as.integer(c(N.start)),
               #= Parameters to save
               beta.nonconst=as.double(beta), ## Fixed parameters of the regression
               gamma.nonconst=as.double(gamma),
               N_pred.nonconst=as.integer(N_pred),
               #= Defining priors
               mubeta=as.double(c(mubeta)), Vbeta=as.double(c(Vbeta)),
               mugamma=as.double(c(mugamma)), Vgamma=as.double(c(Vgamma)),
               #= Diagnostic
               Deviance.nonconst=as.double(Deviance),
               lambda_latent.nonconst=as.double(lambda_latent), ## Predictive posterior mean
               delta_latent.nonconst=as.double(delta_latent), ## Predictive posterior mean
               lambda_pred.nonconst=as.double(lambda_pred),
               #= Seed
               seed=as.integer(seed),
               #= Verbose
               verbose=as.integer(verbose),
               #= Save p and N
               save_p=as.integer(save.p),
               save_N=as.integer(save.N),
               PACKAGE="hSDM")
 
  #= Matrix of MCMC samples
  Matrix <- matrix(NA,nrow=nsamp,ncol=np+nq+1)
  names.fixed <- c(paste("beta.",colnames(X),sep=""),paste("gamma.",colnames(W),sep=""))
  colnames(Matrix) <- c(names.fixed,"Deviance")
  
  #= Filling-in the matrix
  Matrix[,c(1:np)] <- matrix(Sample[[17]],ncol=np)
  Matrix[,c((np+1):(np+nq))] <- matrix(Sample[[18]],ncol=nq)
  Matrix[,ncol(Matrix)] <- Sample[[24]]

  #= Transform Sample list in an MCMC object
  MCMC <- mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)

  #= Save pred
  if (save.p==0) {lambda.pred <- Sample[[27]]}
  if (save.p==1) {
      Matrix.p.pred <- matrix(Sample[[27]],ncol=npred)
      colnames(Matrix.p.pred) <- paste("p.",c(1:npred),sep="")
      lambda.pred <- mcmc(Matrix.p.pred,start=nburn+1,end=ngibbs,thin=nthin)
  }

  #= Save N
  if (save.N==0) {
      N.pred <- Sample[[19]]
  }
  if (save.N==1) {
      Matrix.N.pred <- matrix(Sample[[19]],ncol=nsite)
      colnames(Matrix.N.pred) <- paste("N.",Levels.site,sep="")
      N.pred=mcmc(Matrix.N.pred,start=nburn+1,end=ngibbs,thin=nthin)
  }

  #= Output
  return (list(mcmc=MCMC,
               lambda.pred=lambda.pred, N.pred=N.pred,
               lambda.latent=Sample[[25]], delta.latent=Sample[[26]]))

}

#===================================================================
# END
#===================================================================
