####################################################################
##
## hSDM.siteocc.iCAR.R
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


hSDM.siteocc.iCAR <- function (# Observations
                               presence, observability,
                               site, data.observability,
                               # Habitat
                               suitability, data.suitability,
                               # Spatial structure
                               spatial.entity,
                               n.neighbors, neighbors,
                               # Predictions
                               suitability.pred=NULL, spatial.entity.pred=NULL,
                               # Chains
                               burnin=1000, mcmc=1000, thin=1,
                               # Starting values
                               beta.start,
                               gamma.start,
                               Vrho.start,
                               # Priors
                               mubeta=0, Vbeta=1.0E6,
                               mugamma=0, Vgamma=1.0E6,
                               priorVrho="1/Gamma",
                               shape=0.5, rate=0.0005,
                               Vrho.max=1000,
                               # Various
                               seed=1234, verbose=1,
                               save.rho=0, save.p=0)

{   
  #========
  # Basic checks
  #========
  check.mcmc.parameters(burnin, mcmc, thin)
  check.verbose(verbose)
  check.save.rho(save.rho)
  check.save.p(save.p)
   
  #======== 
  # Form response, covariate matrices and model parameters
  #========

  #= Response
  Y <- presence
  nobs <- length(Y)
  #= Suitability
  mf.suit <- model.frame(formula=suitability,data=data.suitability)
  X <- model.matrix(attr(mf.suit,"terms"),data=mf.suit)
  #= Observability
  mf.obs <- model.frame(formula=observability,data=data.observability)
  W <- model.matrix(attr(mf.obs,"terms"),data=mf.obs)
  #= Sites
  Levels.site <- sort(unique(site))
  nsite <- length(Levels.site)
  sites <- as.numeric(as.factor(site))
  #= Spatial entity
  ncell <- length(n.neighbors)
  cells <- spatial.entity

  #= Predictions
  if (is.null(suitability.pred) | is.null(spatial.entity.pred)) {
      X.pred <- X
      cells.pred <- spatial.entity
      npred <- nsite
  }
  if (!is.null(suitability.pred) & !is.null(spatial.entity.pred)) {
      mf.pred <- model.frame(formula=suitability,data=suitability.pred)
      X.pred <- model.matrix(attr(mf.pred,"terms"),data=mf.pred)
      cells.pred <- spatial.entity.pred
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
  check.Y.binomial(Y,rep(1:nobs))
  check.X(X,nsite) # X must be of dim (nsite x np) for the siteocc.iCAR model
  check.W(W,nobs)
  check.sites(sites,nobs)
  check.cells(cells,nsite)
  check.neighbors(n.neighbors,ncell,neighbors)

  #========
  # Initial starting values for M-H
  #========
  beta.start <- form.beta.start(beta.start,np)
  gamma.start <- form.gamma.start(gamma.start,nq)
  rho.start <- rep(0,ncell) # Starting values for spatial random effects set to zero.
  Vrho.start <- check.Vrho.start(Vrho.start)
  
  #========
  # Form and check priors
  #========
  mubeta <- check.mubeta(mubeta,np)
  Vbeta <- check.Vbeta(Vbeta,np)
  mugamma <- check.mugamma(mugamma,nq)
  Vgamma <- check.Vgamma(Vgamma,nq)
  check.ig.prior(shape,rate)
  Vrho.max <- check.Vrho.max(Vrho.max)
  priorVrho <- form.priorVrho(priorVrho)

  #========
  # Parameters to save
  #========
  beta <- rep(beta.start,nsamp)
  gamma <- rep(gamma.start,nsamp)
  if (save.rho==0) {rho_pred <- rho.start}
  if (save.rho==1) {rho_pred <- rep(rho.start,nsamp)}
  Vrho <- rep(Vrho.start,nsamp)
  theta_latent <- rep(0,nsite)
  delta_latent <- rep(0,nobs)
  if (save.p==0) {theta_pred <- rep(0,npred)}
  if (save.p==1) {theta_pred <- rep(0,npred*nsamp)}
  Deviance <- rep(0,nsamp)

  #========
  # call C++ code to draw sample
  #========
  Sample <- .C("hSDM_siteocc_iCAR",
               #= Constants and data
               ngibbs=as.integer(ngibbs), nthin=as.integer(nthin), nburn=as.integer(nburn), ## Number of iterations, burning and samples
               nobs=as.integer(nobs), nsite=as.integer(nsite), ncell=as.integer(ncell),
               np=as.integer(np),
               nq=as.integer(nq),
               Y_vect=as.integer(c(Y)),
               W_vect=as.double(c(W)),
               X_vect=as.double(c(X)),
               #= Sites
               S_vect=as.integer(c(sites)-1),
               #= Spatial correlation
               C_vect=as.integer(c(cells)-1), # Cells range is 1,...,ncell in R. Must start at 0 for C. Don't forget the "-1" term. 
               nNeigh=as.integer(c(n.neighbors)),
               Neigh_vect=as.integer(c(neighbors-1)), # Cells range is 1,...,ncell in R. Must start at 0 for C. Don't forget the "-1" term.
               #= Predictions
               npred=as.integer(npred),
               X_pred_vect=as.double(c(X.pred)),
               C_pred_vect=as.integer(c(cells.pred)-1),
               #= Starting values for M-H
               beta_start=as.double(c(beta.start)),
               gamma_start=as.double(c(gamma.start)),
               rho_start=as.double(c(rho.start)),
               #= Parameters to save
               beta.nonconst=as.double(beta), ## Fixed parameters of the regression
               gamma.nonconst=as.double(gamma),
               rho_pred.nonconst=as.double(rho_pred), 
               Vrho.nonconst=as.double(Vrho),
               #= Defining priors
               mubeta=as.double(c(mubeta)), Vbeta=as.double(c(Vbeta)),
               mugamma=as.double(c(mugamma)), Vgamma=as.double(c(Vgamma)),
               priorVrho=as.double(priorVrho),
               shape=as.double(shape), rate=as.double(rate),
               Vrho.max=as.double(Vrho.max),
               #= Diagnostic
               Deviance.nonconst=as.double(Deviance),
               theta_latent.nonconst=as.double(theta_latent), ## Predictive posterior mean
               delta_latent.nonconst=as.double(delta_latent), ## Predictive posterior mean
               theta_pred.nonconst=as.double(theta_pred),
               #= Seed
               seed=as.integer(seed),
               #= Verbose
               verbose=as.integer(verbose),
               #= Save rho and p
               save_rho=as.integer(save.rho),
               save_p=as.integer(save.p),
               PACKAGE="hSDM")
 
  #= Matrix of MCMC samples
  Matrix <- matrix(NA,nrow=nsamp,ncol=np+nq+2)
  names.fixed <- c(paste("beta.",colnames(X),sep=""),paste("gamma.",colnames(W),sep=""))
  colnames(Matrix) <- c(names.fixed,"Vrho","Deviance")
  
  #= Filling-in the matrix
  Matrix[,c(1:np)] <- matrix(Sample[[22]],ncol=np)
  Matrix[,c((np+1):(np+nq))] <- matrix(Sample[[23]],ncol=nq)
  Matrix[,ncol(Matrix)-1] <- Sample[[25]]
  Matrix[,ncol(Matrix)] <- Sample[[34]]

  #= Transform Sample list in an MCMC object
  MCMC <- mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)

  #= Save rho
  if (save.rho==0) {rho.pred <- Sample[[24]]}
  if (save.rho==1) {
      Matrix.rho.pred <- matrix(Sample[[24]],ncol=ncell)
      colnames(Matrix.rho.pred) <- paste("rho.",c(1:ncell),sep="")
      rho.pred <- mcmc(Matrix.rho.pred,start=nburn+1,end=ngibbs,thin=nthin)
  }
  
  #= Save pred
  if (save.p==0) {theta.pred <- Sample[[37]]}
  if (save.p==1) {
      Matrix.theta.pred <- matrix(Sample[[37]],ncol=npred)
      colnames(Matrix.theta.pred) <- paste("p.",c(1:npred),sep="")
      theta.pred <- mcmc(Matrix.theta.pred,start=nburn+1,end=ngibbs,thin=nthin)
  }

  #= Output
  return (list(mcmc=MCMC,
               rho.pred=rho.pred, theta.pred=theta.pred,
               theta.latent=Sample[[35]], delta.latent=Sample[[36]]))

}

#===================================================================
# END
#===================================================================
