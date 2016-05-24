# ***
# Functions to extract LR statistic for supported classes of statistical models,
# possibly given user-defined values for model parameters
#

# *** lmer method
.getLR.lmer <- function(model, null.model=NULL, psi=NULL, null.psi=NULL){

  if(!requireNamespace("lme4", quietly=TRUE)) stop("The 'lme4' package must be installed in order to use this function.")
  requireNamespace("lme4", quietly=TRUE)

  if(is.null(null.model)){
  # TODO: global model fit (experimental, for D2)

    k <- attr(logLik(model[[1]]),"df")
    reml <- lme4::isREML(model[[1]])
    logL1 <- sapply(model,logLik)
    dW <- -2*logL1

  }else{

    k <- attr(logLik(model[[1]]),"df") - attr(logLik(null.model[[1]]),"df") 
    reml <- any(lme4::isREML(model[[1]]), lme4::isREML(null.model[[1]]))
    if(is.null(psi) & is.null(null.psi)){
      logL0 <- sapply(null.model,logLik)
      logL1 <- sapply(model,logLik)
    }else{
      logL0 <- sapply(null.model,.logL.lmer,psi=null.psi)
      logL1 <- sapply(model,.logL.lmer,psi=psi)
    }
    dW <- -2*(logL0-logL1)

  }
  attr(dW,"df") <- k
  attr(dW,"REML") <- reml
  dW

}

# *** nlme method
.getLR.nlme <- function(model, null.model=NULL, psi=NULL, null.psi=NULL){

  if(!requireNamespace("nlme", quietly=TRUE)) stop("The 'nlme' package must be installed in order to use this function.")
  requireNamespace("nlme", quietly=TRUE)

  if(is.null(null.model)){
  # TODO: global model fit (experimental, for D2)

    k <- attr(logLik(model[[1]]),"df")
    reml <- model[[1]]$method=="REML"
    logL1 <- sapply(model,logLik)
    dW <- -2*logL1

  }else{

    k <- attr(logLik(model[[1]]),"df") - attr(logLik(null.model[[1]]),"df") 
    reml <- any(model[[1]]$method=="REML", null.model[[1]]$method=="REML")
    if(is.null(psi) & is.null(null.psi)){
      logL0 <- sapply(null.model,logLik)
      logL1 <- sapply(model,logLik)
    }else{
      logL0 <- sapply(null.model,.logL.nlme,psi=null.psi)
      logL1 <- sapply(model,.logL.nlme,psi=psi)
    }
    dW <- -2*(logL0-logL1)

  }
  attr(dW,"df") <- k
  attr(dW,"REML") <- reml
  dW

}

# *** default method
.getLR.default <- function(model,null.model=NULL, psi=NULL, null.psi=NULL){


  if(is.null(null.model)){
  # TODO: global model fit (experimental, D2)

    k <- .tryResidualDf(model[[1]])
    logL1 <- sapply(model,logLik)
    dW <- -2*logL1

  }else{

    k <- .tryResidualDf(null.model[[1]]) - .tryResidualDf(model[[1]])
    if(is.null(psi) & is.null(null.psi)){
      logL0 <- sapply(null.model,logLik)
      logL1 <- sapply(model,logLik)
    }else{
      # only lm supported
      logL0 <- sapply(null.model,.logL.lm,psi=null.psi)
      logL1 <- sapply(model,.logL.lm,psi=psi)
    }
    dW <- -2*(logL0-logL1)

  }
  attr(dW,"df") <- k
  dW

}


# ***
# Likelihood functions for supported models
#

.logL.lmer <- function(object, psi=NULL){

  if(is.null(psi)){
    beta <- lme4::fixef(object)
    D <- lme4::VarCorr(object)[[1]]
    sig2 <- attr(lme4::VarCorr(object),"sc")^2
  }else{
    beta <- psi$beta
    D <- psi$D
    sig2 <- psi$sigma2
  }

  cls <- lme4::getME(object,"flist")[[1]]
  p <- length(beta)
  q <- dim(D)[1]
  y <- split(lme4::getME(object,"y"),cls)
  X <- split(lme4::getME(object,"X"),cls)
  # testing: via mmList, for single level of clustering
  Z <- split(lme4::getME(object,"mmList")[[1]],cls)
  # outdated: via sparse model matrix (produces bug due to zero elements)
  # Z <- split(matrix(lme4::getME(object,"Zt")@x, nrow=length(cls), byrow=T),cls)

  L <- numeric(nlevels(cls))
  for(i in levels(cls)){

    yi <- y[[i]]
    ni <- length(yi)
    Xi <- matrix(X[[i]], ncol=p)
    Ri <- yi - Xi%*%beta
    Zi <- matrix(Z[[i]], ncol=q)

    V <- diag(sig2,ni) + Zi %*% D %*% t(Zi)
    Vinv <- chol2inv(chol(V))
    
    dV <- determinant(V,logarithm=TRUE)
    dV <- dV$modulus*dV$sign
    L[i] <- dV + t(Ri) %*% Vinv %*% (Ri)
  }
  -sum(L)/2

}

.logL.nlme <- function(object, psi=NULL){

  if(is.null(psi)){
    beta <- nlme::fixef(object)
    D <- nlme::getVarCov(object)
    sig2 <- object$sigma^2
  }else{
    beta <- psi$beta
    D <- psi$D
    sig2 <- psi$sigma2
  }

  # error check
  if(is.null(nlme::getData(object))) stop("No data sets found in 'lme' fit. See '?testModels' for an example.")

  cls <- nlme::getGroups(object)
  p <- length(beta)
  q <- dim(D)[1]
  y <- split(nlme::getResponse(object),cls)
  fe <- object$terms
  X <- split(model.matrix(fe,nlme::getData(object)),cls)
  re <- attr(object$modelStruct$reStruct[[1]],"formula")
  Z <- split(model.matrix(re,nlme::getData(object)),cls)

  L <- numeric(nlevels(cls))
  for(i in levels(cls)){

    yi <- y[[i]]
    ni <- length(yi)
    Xi <- matrix(X[[i]], ncol=p)
    Ri <- yi - Xi%*%beta
    Zi <- matrix(Z[[i]], ncol=q)

    V <- diag(sig2,ni) + Zi %*% D %*% t(Zi)
    Vinv <- chol2inv(chol(V))
    
    dV <- determinant(V,logarithm=TRUE)
    dV <- dV$modulus*dV$sign
    L[i] <- dV + t(Ri) %*% Vinv %*% (Ri)
  }
  -sum(L)/2

}

.logL.lm <- function(object, psi=NULL){

  if(is.null(psi)){
    beta <- coef(object)
    # sig2 <- sum(resid(object)^2)/df.residual(object)
    r <- resid(object)                # SiG 19-04-2016
    sig2 <- sum(r^2)/length(r)
  }else{
    beta <- psi$beta
    sig2 <- psi$sigma2
  }
  ytrm <- attr(object$terms,"variables")[-1][attr(object$terms,"response")]
  y <- as.matrix(object$model[as.character(ytrm)])
  X <- model.matrix(object$terms,object$model)
  r <- y-X%*%beta
  n <- length(y)

  L <- -(n/2)*log(sig2) - (1/(2*sig2)) * t(r) %*% r
  as.numeric(L)
  
}

.tryResidualDf <- function(object){

    k <- NULL
    if(is.null(k)) k <- tryCatch( df.residual(object), error=function(f) NULL )
    if(is.null(k)) k <- tryCatch( tail(anova(object),1), error=function(f) NULL )
    k

}
