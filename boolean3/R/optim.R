## ----------------------------------------------------------------------------
## This file is part of boolean3
##
## Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>
##
## boolean3 represents a substantial re-write of the original boolean package
## developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
## was developed under the direction of Bear Braumoeller and with support from
## The Ohio State University's College of Social and Behavioral Sciences.
##
## boolean3 is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free
## Software Foundation, either version 3 of the License, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program.  If not, see <http://www.gnu.org/licenses/>.
##
## ----------------------------------------------------------------------------

make_lik <- function(obj, positive=FALSE) {
  ## Build the log-likelihood function for latter evaluation. Returns a
  ## closure.
  
  ## Set up the submodel calc. This is not evaluated here, but rather when
  ## Fn calls eval.
  for (mod in names(obj$link)) {
    if (obj$link[[mod]]$name == "scobit") {
      assign (mod,
              substitute (obj$link[[mod]]$linkinv
                          (tcrossprod
                           (beta[obj$coef.idx[[mod]]],
                            obj$frame[[mod]]),
                           shape=beta[paste("lambda", mod, sep="_")])))
    }
    else {
      assign (mod,
              substitute (obj$link[[mod]]$linkinv
                          (tcrossprod
                           (beta[obj$coef.idx[[mod]]],
                            obj$frame[[mod]]))))
    }
  }

  y <- obj$response

  if (positive)
    llik <- substitute(sum(y*(log.p - log.p1) + log.p1))
  else
    llik <- substitute(-sum(y*(log.p - log.p1) + log.p1))

  ## Define closure.
  Fn <- function(beta) {
    for (n in names(obj$model$submods))
      assign(n, eval(get(n)))
    p <- eval(obj$model$boolmod)
    log.p  <- log(p)
    log.p1 <- log(1-p)
    eval(llik)
  }

  return(Fn)
}

make_prob <- function(obj) {
  ## Build the predicted probability function for later evaluation. Returns a
  ## closure.

  beta <- obj$model.fit[[1]]$par

  ## Set up the submodel calc. This is not evaluated here, but rather when
  ## Fn calls eval.
  for (mod in names(obj$link)) {
    if (obj$link[[mod]]$name == "scobit") {
      assign (mod,
              substitute (obj$link[[mod]]$linkinv
                          (crossprod
                           (beta[obj$coef.idx[[mod]]],
                            x[obj$coef.idx[[mod]]]),
                           shape=beta[paste("lambda", mod, sep="_")])))
    }
    else {      
      assign (mod,
              substitute (obj$link[[mod]]$linkinv
                          (crossprod
                           (beta[obj$coef.idx[[mod]]],
                            x[obj$coef.idx[[mod]]]))))
    }
  }
  
  ## Define closure.
  Fn <- function(x) {
    for (n in names(obj$model$submods))
      assign(n, eval(get(n)))
    eval(obj$model$boolmod)
  }

  return(Fn)
}

calc_hessian <- function(mle, obj) {
  ## Calculate hessian using numDeriv's numerical algorithm.
  hessian(obj$calc_lik, mle)
}

set_start <- function(obj, ...) {
  ## Calculate coefficient starting values.
  
  ## Construct the matrix
  F <- matrix(obj$response, ncol=1, dimnames=list(NULL, "y"))
  M <- c()
  for (n in names(obj$frame)) {
    F <- cbind(F, obj$frame[[n]])
    M <- c(M, obj$model$submods[[n]])
  }
  F <- as.data.frame(F)

  ## Fit the model.
  b <- coef(glm(paste("y ~", paste(M, collapse = "+")), data = F))

  ## Build starting values.
  B <- c()
  for (n in names(obj$frame))
    B <- c(B, b[colnames(obj$frame[[n]])]) 

  if (!(is.null(obj$scobit.lambda)))
    B <- c(B, obj$scobit.lambda)

  B
}

fit_optimx <- function(obj, start=start, ...)
{
    ## Fit model with any of the methods available in optimx.
    method <- obj$method
    fit <- optimx(start, obj$calc_lik, hessian=TRUE, method=method, ...)
    details <- attr(fit, "details")
    
    keep <- 1:(ncol(fit) - 8)
    result <- list()

    for (mod in rownames(fit)) {
        est.par <- as.matrix(fit[mod, keep])        
        est.hessian <- details[mod,]$nhatend
        result[[mod]] <- list(par=est.par,
                              nhatend=est.hessian,
                              value=-1*fit[mod,"value"],
                              iterations=fit[mod,"niter"],
                              conv=fit[mod,"convcode"],
                              detail=fit)        
    }

    result
}

fit_optim <- function(obj, start=start, ...) {
    ## Fit model with simulated-annealing, a method available in optim but not
    ## in optimx.
    fit <- optim(start, obj$calc_lik, method="SANN", hessian=FALSE, ...)
    fit$hessian <- calc_hessian(fit$par, obj)

    list(SANN=list(par=fit$par, nhatend=fit$hessian, value=-1*fit$value,
             iterations=fit$counts[1], conv=0, SANN=fit))
}

fit_genoud <- function(obj, start=start, ...) {
    ## Fit model with the genetic optimization algorithm available in genoud.
  
    ## If a scobit model has been specified, boundary constraints need to be set
    ## up to guarantee lambda > 0.
    Domains <- list(...)$Domains
    if (!is.null(obj$scobit.lambda)) {
        Domains <- matrix(0, ncol=2, nrow=obj$k)
    }
  
    fit <- genoud(obj$calc_lik, obj$k, starting.values=start, hessian=TRUE,
                  Domains=Domains, ...)

    list(genoud=list(par=fit$par, nhatend=fit$hessian, value=-1*fit$value,
             iterations=fit$generations, conv=0, genoud=fit))
}

## fit_mcmc <- function(obj, start=start, ...) {
##   ## Fit model with MCMC from MCMCpack.

##   ## Get initial fit with optimx.
##   tmp <- obj
##   tmp$method <- "nlminb"
##   initial_fit <- fit_optimx(tmp, start=start)

##   ## likelihood needs to be positive for this function.
##   obj$calc_lik <- make_lik(obj, positive=TRUE)
##   fit <- MCMCmetrop1R(obj$calc_lik, initial_fit$nlminb$par, verbose=500, 
##                       V=solve(initial_fit$nlminb$nhatend), ...)
                      
##   smry <- summary(fit)

##   list(mcmc=list(coef=smry$statistics[,1], chains=smry$nchain, thin=smry$thin,
##        iterations=smry$end - smry$start + 1, mcmc=fit))
## }
