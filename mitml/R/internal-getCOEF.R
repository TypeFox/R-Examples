# ***
# Functions to extract (fixed) coefficients and SEs/Covariance from 
# supported classes of statistical models
#

# *** lmer method
.getCOEF.lmer <- function(model,null.model=NULL,diagonal=FALSE){

  if(!requireNamespace("lme4", quietly=TRUE)) stop("The 'lme4' package must be installed in order to use this function.")
  if(is.null(null.model)){

    if(diagonal){
      nms <- names(lme4::fixef(model[[1]]))
      Qhat <- sapply(model,lme4::fixef)
      if(is.null(dim(Qhat))){ Uhat <- sapply(model, function(z) vcov(summary(z))@x )
      }else{ Uhat <- sapply(model, function(z) diag( matrix(vcov(summary(z))@x,nrow(Qhat)) )) }
    }else{
      nms <- names(lme4::fixef(model[[1]]))
      p <- length(nms)
      Qhat <- sapply(model,lme4::fixef)
      Uhat <- vapply(model, function(z) vcov(summary(z))@x, FUN.VALUE=matrix(0,p,p))
    }

  }else{

    par0 <- names(lme4::fixef(null.model[[1]]))
    par1 <- names(lme4::fixef(model[[1]]))
    dpar <- setdiff(par1,par0)
    nms <- NULL
    p <- length(par1)
    i <- which(par1%in%dpar)
    Qhat <- sapply(model,lme4::fixef)[dpar,]
    Uhat <- vapply(model, function(z) vcov(summary(z))@x, FUN.VALUE=matrix(0,p,p))[i,i,]

  }
  out <- list(Qhat=Qhat,Uhat=Uhat,nms=nms)
  out

}

# *** nlme method
.getCOEF.nlme <- function(model,null.model=NULL,diagonal=FALSE){

  if(!requireNamespace("nlme", quietly=TRUE)) stop("The 'nlme' package must be installed in order to use this function.")
  if(is.null(null.model)){

    if(diagonal){
      nms <- names(nlme::fixef(model[[1]]))
      Qhat <- sapply(model,nlme::fixef)
      if(is.null(dim(Qhat))){ Uhat <- sapply(model, function(z) vcov(summary(z)) )
      }else{ Uhat <- sapply(model, function(z) diag( vcov(summary(z)) ) ) }
    }else{
      nms <- names(nlme::fixef(model[[1]]))
      p <- length(nms)
      Qhat <- sapply(model,nlme::fixef)
      Uhat <- vapply(model, function(z) vcov(summary(z)), FUN.VALUE=matrix(0,p,p))
    }

  }else{

    par0 <- names(nlme::fixef(null.model[[1]]))
    par1 <- names(nlme::fixef(model[[1]]))
    dpar <- setdiff(par1,par0)
    nms <- NULL
    p <- length(par1)
    i <- which(par1%in%dpar)
    Qhat <- sapply(model,nlme::fixef)[dpar,]
    Uhat <- vapply(model, function(z) vcov(summary(z)), FUN.VALUE=matrix(0,p,p))[i,i,]

  }
  out <- list(Qhat=Qhat,Uhat=Uhat,nms=nms)
  out

}

# *** default method
.getCOEF.default <- function(model,null.model=NULL,diagonal=FALSE){

  if(is.null(null.model)){

    if(diagonal){
      nms <- names(coef(model[[1]]))
      Qhat <- sapply(model,coef)
      if(is.null(dim(Qhat))){ Uhat <- sapply(model, vcov )
      }else{ Uhat <- sapply(model,function(z) diag(vcov(z))) }
    }else{
      nms <- names(coef(model[[1]]))
      p <- length(nms)
      Qhat <- sapply(model,coef)
      Uhat <- vapply(model,vcov, FUN.VALUE=matrix(0,p,p))
    }

  }else{

    par0 <- names(coef(null.model[[1]]))
    par1 <- names(coef(model[[1]]))
    dpar <- setdiff(par1,par0)
    nms <- NULL
    p <- length(par1)
    i <- which(par1%in%dpar)
    Qhat <- sapply(model,coef)[dpar,]
    Uhat <- vapply(model,vcov, FUN.VALUE=matrix(0,p,p))[i,i,]
 
  }
  out <- list(Qhat=Qhat,Uhat=Uhat,nms=nms)
  out

}

