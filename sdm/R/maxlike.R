# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
.maxlike <- function(formula,data,link='logit',starts,hessian=TRUE,fixed,normalize=TRUE,...) {
  # adopted based on the function maxlike in the maxlike package by Richard Chandler and Andy Royle
  
  link <- .pmatch(link,c("logit","cloglog"))[1]
  if (is.na(link)) stop("link function should be either 'logit' or 'cloglog'")
  varnames <- all.vars(formula)
  call <- match.call()
  
  nsp <- deparse(formula[[2]])
  if (normalize) data <- .normalize(data,except=nsp)
  w <- data[,nsp] == 1
  X <- data[w,]
  Z <- data[-w,]
  
  npts <- nrow(X)
  X <- model.frame(formula,X,na.action = 'na.omit')
  X <- model.matrix(formula, X)
  
  Z <- model.frame(formula,Z,na.action = 'na.omit')
  Z <- model.matrix(formula, Z)
  
  
  npars <- ncol(X)
  parnames <- colnames(X)
  
  if (missing(starts) || length(starts) != length(npars)) {
    starts <- rep(0, npars)
    names(starts) <- parnames
  } else names(starts) <- parnames
  
  if (identical(link, "logit")) {
    nll <- function(pars) {
      psix <- plogis(drop(X %*% pars))
      psiz <- sum(plogis(drop(Z %*% pars)))
      -1 * sum(log(psix/psiz))
    }
  } else if (identical(link, "cloglog")) {
    nll <- function(pars) {
      psix <- 1 - exp(-exp(drop(X %*% pars)))
      psiz <- sum(1 - exp(-exp(drop(Z %*% pars))))
      -1 * sum(log(psix/psiz))
    }
  } else stop("link function should be either 'logit' or 'cloglog'")
  
  is.fixed <- rep(FALSE, npars)
  
  if (!missing(fixed)) {
    if (length(fixed) != length(starts))
      stop("fixed should be a vector with the same length as the number of parameters to be estimated")
    if (sum(is.double(fixed)) < 1) 
      stop("fixed must contain at least one real value")
    is.fixed <- !is.na(fixed)
    if (sum(!is.fixed) < 1) 
      stop("you cannot fix all parameters in the model")
    npars <- sum(!is.fixed)
    nll.fix <- function(p) {
      p[is.fixed] <- fixed[is.fixed]
      do.call("nll", list(pars = p))
    }
    fm <- optim(starts, nll.fix, hessian = hessian, ...)
    fm$par[is.fixed] <- fixed[is.fixed]
  } else {
    fm <- optim(starts, nll, hessian = hessian, ...)
  }
  not.fixed <- !is.fixed
  par <- fm$par
  if (hessian) {
    vcTry <- try(solve(fm$hessian[not.fixed, not.fixed]),silent=TRUE)
    if (identical(class(vcTry), "matrix")) {
      vc <- matrix(0, length(par), length(par))
      vc[not.fixed, not.fixed] <- vcTry
      se <- sqrt(diag(vc))
    }
    else {
      vc <- matrix(NA, npars, npars)
      se <- rep(NA, npars)
    }
  }
  else {
    vc <- matrix(NA, npars, npars)
    se <- rep(NA, npars)
  }
  dimnames(vc) <- list(parnames, parnames)
  aic <- 2 * fm$value + 2 * npars
  new('.maxlikeModel',fit=list(Est = cbind(Est = par, SE = se), vcov = vc, 
                               AIC = aic, call = call, optim = fm, not.fixed = not.fixed, link = link,normalize=normalize))
  
}


# adopted based on the function predict.maxlikeFit in the maxlike package by Richard Chandler and Andy Royle
setMethod('predict', signature(object='.maxlikeModel'), 
          function(object, newdata,...) {
            if (missing(newdata)) stop('mewdata is missing...')
            e <- object@fit$Est[, "Est"]
            
            link <- object@fit$link
            npix <- nrow(newdata)
            formula <- .getRhsFromFormula(object@fit$call$formula)
            varnames <- all.vars(formula)
            if (!all(varnames %in% colnames(newdata))) stop("at least 1 covariate in the formula is not in rasters.")
            if (object@fit$normalize) newdata <- .normalize(newdata[,varnames])
            Z <- model.frame(formula, newdata, na.action = "na.pass")
            Z.terms <- attr(Z, "terms")
            Z <- model.matrix(Z.terms, Z)
            
            eta <- drop(Z %*% e)
            if (identical(link, "logit")) 
              eta <- plogis(eta)
            else if (identical(link, "cloglog")) 
              eta <- 1 - exp(-exp(eta))
            else stop("link function should be either 'logit' or 'cloglog'")
            eta
          }
)
