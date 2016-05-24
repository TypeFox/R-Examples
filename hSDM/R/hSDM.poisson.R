###################################################################
##
## hSDM.poisson.R
##
####################################################################
##
## Original code by Ghislain Vieilledent, October 2011
## CIRAD UR B&SEF
## ghislain.vieilledent@cirad.fr / ghislainv@gmail.com
##
####################################################################
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991. See the package LICENSE
## file for more information.
##
## Copyright (C) 2011 Ghislain Vieilledent
## 
####################################################################


hSDM.poisson <- function (counts,
                          suitability,
                          data,
                          suitability.pred=NULL,
                          burnin=5000, mcmc=10000, thin=10, 
                          beta.start,
                          mubeta=0, Vbeta=1.0E6,
                          seed=1234, verbose=1, save.p=0)

{   
  #========
  # Basic checks
  #========
  check.mcmc.parameters(burnin, mcmc, thin)
  check.verbose(verbose)
  check.save.p(save.p)
   
  #======== 
  # Form response, covariate matrices and model parameters
  #========

  #= Response
  Y <- counts
  nobs <- length(Y)
  #= Suitability
  mf.suit <- model.frame(formula=suitability,data=data)
  X <- model.matrix(attr(mf.suit,"terms"),data=mf.suit)
  #= Predictions
  if (is.null(suitability.pred)) {
      X.pred <- X
      npred <- nobs
  }
  if (!is.null(suitability.pred)) {
      mf.pred <- model.frame(formula=suitability,data=suitability.pred)
      X.pred <- model.matrix(attr(mf.pred,"terms"),data=mf.pred)
      npred <- nrow(X.pred)
  }
  #= Model parameters
  np <- ncol(X)
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin

  #========== 
  # Check data
  #==========
  check.Y.poisson(Y)
  check.X(X,nobs)
  
  #========
  # Initial starting values for M-H
  #========
  beta.start <- form.beta.start(beta.start,np)

  #========
  # Form and check priors
  #========
  mubeta <- check.mubeta(mubeta,np)
  Vbeta <- check.Vbeta(Vbeta,np)

  #========
  # Parameters to save
  #========
  beta <- rep(beta.start,nsamp)
  prob_p_latent <- rep(0,nobs)
  if (save.p==0) {prob_p_pred <- rep(0,npred)}
  if (save.p==1) {prob_p_pred <- rep(0,npred*nsamp)}
  Deviance <- rep(0,nsamp)

  #========
  # call C++ code to draw sample
  #========
  Sample <- .C("hSDM_poisson",
               #= Constants and data
               ngibbs=as.integer(ngibbs), nthin=as.integer(nthin), nburn=as.integer(nburn), ## Number of iterations, burning and samples
               nobs=as.integer(nobs),
               np=as.integer(np),
               Y_vect=as.integer(c(Y)),
               X_vect=as.double(c(X)),
               #= Predictions
               npred=as.integer(npred),
               X_pred_vect=as.double(c(X.pred)),
               #= Starting values for M-H
               beta_start=as.double(c(beta.start)),
               #= Parameters to save
               beta.nonconst=as.double(beta), ## Fixed parameters of the regression
               #= Defining priors
               mubeta=as.double(c(mubeta)), Vbeta=as.double(c(Vbeta)),
               #= Diagnostic
               Deviance.nonconst=as.double(Deviance),
               prob_p_latent.nonconst=as.double(prob_p_latent), ## Predictive posterior mean
               prob_p_pred.nonconst=as.double(prob_p_pred), 
               #= Seed
               seed=as.integer(seed), 
               #= Verbose
               verbose=as.integer(verbose),
               #= Save p
               save_p=as.integer(save.p),
               PACKAGE="hSDM")
 
  #= Matrix of MCMC samples
  Matrix <- matrix(NA,nrow=nsamp,ncol=np+1)
  names.fixed <- paste("beta.",colnames(X),sep="")
  colnames(Matrix) <- c(names.fixed,"Deviance")
  
  #= Filling-in the matrix
  Matrix[,c(1:np)] <- matrix(Sample[[11]],ncol=np)
  Matrix[,ncol(Matrix)] <- Sample[[14]]

  #= Transform Sample list in an MCMC object
  MCMC <- mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)

  #= Save pred
  if (save.p==0) {lambda.pred <- Sample[[16]]}
  if (save.p==1) {
      Matrix.p.pred <- matrix(Sample[[16]],ncol=npred)
      colnames(Matrix.p.pred) <- paste("p.",c(1:npred),sep="")
      lambda.pred <- mcmc(Matrix.p.pred,start=nburn+1,end=ngibbs,thin=nthin)
  }
  
  #= Output
  return (list(mcmc=MCMC, lambda.pred=lambda.pred, lambda.latent=Sample[[15]]))

}

#===================================================================
# END
#===================================================================
