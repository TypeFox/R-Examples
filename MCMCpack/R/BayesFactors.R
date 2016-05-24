##########################################################################
## BayesFactor.R contains functions useful for calculating and comparing
## marginal likelihoods
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Originally written by KQ 1/26/2006
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


## log densities
"logdinvgamma" <- function(sigma2, a, b){
  logf <- a * log(b) - lgamma(a) + -(a+1) * log(sigma2) + -b/sigma2
  return(logf)
}

"logdmvnorm" <- function(theta, mu, Sigma){
  d <- length(theta)
  logf <- -0.5*d * log(2*pi) - 0.5*log(det(Sigma)) -
    0.5 * t(theta - mu) %*% solve(Sigma) %*% (theta - mu)
  return(logf)
}






## log posterior densities
"logpost.regress" <- function(theta, y, X, b0, B0, c0, d0){
  n <- length(y)
  k <- ncol(X)
  beta <- theta[1:k]
  sigma2 <- exp(theta[k+1])

  Sigma <- solve(B0)
  
  loglike <- sum(dnorm(y, X%*%beta, sqrt(sigma2), log=TRUE))

  ## note the change to the prior for sigma2 b/c of the transformation
  logprior <- logdinvgamma(sigma2, c0/2, d0/2) + theta[k+1] +
    logdmvnorm(beta, b0, Sigma)

  return (loglike + logprior)  
}


"logpost.probit" <- function(theta, y, X, b0, B0){
  n <- length(y)
  k <- ncol(X)
  beta <- theta

  p <- pnorm(X %*% beta)
  
  Sigma <- solve(B0)
  
  loglike <- sum( y * log(p) + (1-y)*log(1-p) )
  logprior <- logdmvnorm(beta, b0, Sigma)

  return (loglike + logprior)  
}

"logpost.logit" <- function(theta, y, X, b0, B0){
  n <- length(y)
  k <- ncol(X)
  beta <- theta

  eta <- X %*% beta
  p <- 1.0/(1.0+exp(-eta))
  
  Sigma <- solve(B0)
  
  loglike <- sum( y * log(p) + (1-y)*log(1-p) )
  logprior <- logdmvnorm(beta, b0, Sigma)

  return (loglike + logprior)  
}

"logpost.logit.userprior" <- function(theta, y, X, userfun, logfun,
                                      my.env){
  n <- length(y)
  k <- ncol(X)
  beta <- theta

  eta <- X %*% beta
  p <- 1.0/(1.0+exp(-eta))
  
  loglike <- sum( y * log(p) + (1-y)*log(1-p) )
  if (logfun){
    logprior <- eval(userfun(theta), envir=my.env)
  }
  else{
    logprior <- log(eval(userfun(theta), envir=my.env))
  }

  return (loglike + logprior)  
}


"logpost.poisson" <- function(theta, y, X, b0, B0){
  n <- length(y)
  k <- ncol(X)
  beta <- theta

  eta <- X %*% beta
  lambda <- exp(eta)
  
  Sigma <- solve(B0)
  
  loglike <- sum(dpois(y, lambda, log=TRUE))
  logprior <- logdmvnorm(beta, b0, Sigma)
  
  return (loglike + logprior)  
}




## functions for working with BayesFactor objects
"BayesFactor" <- function(...){
  model.list <- list(...)
  M <- length(model.list)  
  #model.names <- paste("model", 1:M, sep="")
  this.call <- match.call()
  this.call.string <- deparse(this.call)
  this.call.string <- strsplit(this.call.string, "BayesFactor\\(")
  this.call.string <- this.call.string[[1]][length(this.call.string[[1]])]
  this.call.string <- strsplit(this.call.string, ",")

  model.names <- NULL
  for (i in 1:M){
    model.names <- c(model.names, this.call.string[[1]][i])
  }
  model.names <- gsub(")", "", model.names)
  model.names <- gsub(" ", "", model.names)
  
  for (i in 1:M){
    if (!is.mcmc(model.list[[i]])){
      stop("argument not of class mcmc\n")
    }
  }
  
  BF.mat <- matrix(NA, M, M)
  BF.log.mat <- matrix(NA, M, M)
  rownames(BF.mat) <- colnames(BF.mat) <-
      rownames(BF.log.mat) <- colnames(BF.log.mat) <- model.names

  BF.logmarglike <- array(NA, c(1, M), dimnames=list("logmarglike", model.names))
  ## Bill Dunlap found that R did not warn about illegal dimnames in array()  but rather would just disregard them.
  ## So based on the patch by Martin Maechler, JHP changed the code. "Thu Feb  4 10:35:15 2016"
  ## BF.logmarglike <- array(NA, M, dimnames=model.names)
  BF.call <- NULL
  
  for (i in 1:M){
    BF.logmarglike[i] <- attr(model.list[[i]], "logmarglike")
    BF.call <- c(BF.call, attr(model.list[[i]], "call"))
    for (j in 1:M){
      if (identical(attr(model.list[[i]], "y"), attr(model.list[[j]], "y"))){
        BF.log.mat[i,j] <- attr(model.list[[i]], "logmarglike") -
          attr(model.list[[j]], "logmarglike")
        BF.mat[i,j] <- exp(BF.log.mat[i,j])
      }
      else{
          warning(paste(model.names[i], " and ", model.names[j],
                        " do not have exactly identical y data.\nBayes factors are not defined.\n", sep=""))
      }
    }
  }
  
  return(structure(list(BF.mat=BF.mat, BF.log.mat=BF.log.mat,
                        BF.logmarglike=BF.logmarglike,
                        BF.call=BF.call),
                   class="BayesFactor"))
}


"is.BayesFactor" <- function(BF){
  return(class(BF) == "BayesFactor")
}


"print.BayesFactor" <- function(x, ...){

  cat("The matrix of Bayes Factors is:\n")
  print(x$BF.mat, digits=3)

  cat("\nThe matrix of the natural log Bayes Factors is:\n")
  print(x$BF.log.mat, digits=3)

  M <- length(x$BF.call)
  for (i in 1:M){
    cat("\n", rownames(x$BF.mat)[i], ":\n")
    cat("   call = \n")
    print(x$BF.call[[i]])
    cat("\n   log marginal likelihood = ", x$BF.logmarglike[i], "\n\n")
    
  }
  
}


"summary.BayesFactor" <- function(object, ...){

  cat("The matrix of Bayes Factors is:\n")
  print(object$BF.mat, digits=3)

  cat("\nThe matrix of the natural log Bayes Factors is:\n")
  print(object$BF.log.mat, digits=3)

  BF.mat.NA <- object$BF.mat    
  diag(BF.mat.NA) <- NA
  minvec <- apply(BF.mat.NA, 1, min, na.rm=TRUE)
  best.model <- which.max(minvec)
  if (minvec[best.model] > 150){
    cat("\nThere is very strong evidence to support",
        rownames(object$BF.mat)[best.model],
        "over\nall other models considered.\n") 
  }
  else if(minvec[best.model] > 20){
    cat("\nThere is strong evidence or better to support",
        rownames(object$BF.mat)[best.model],
        "over\nall other models considered.\n")    
  }
  else if(minvec[best.model] > 3){
    cat("\nThere is positive evidence or better to support",
        rownames(object$BF.mat)[best.model],
        "over\nall other models considered.\n")    
  }
  else {
    cat("\nThe evidence to support",
        rownames(object$BF.mat)[best.model],
        "over all\nother models considered is worth no more\n than a bare mention.\n")    
  }

  cat("\n\nStrength of Evidence Guidelines\n(from Kass and Raftery, 1995, JASA)\n")
  cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
  cat("2log(BF[i,j])       BF[i,j]         Evidence Against Model j\n")
  cat("------------------------------------------------------------\n")
  cat("  0 to 2            1 to 3           Not worth more than a\n")
  cat("                                          bare mention\n")
  cat("  2 to 6            3 to 20          Positive\n")
  cat("  6 to 10           20 to 150        Strong\n")
  cat("  >10               >150             Very Strong\n")
  cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n")
  
  M <- length(object$BF.call)
  for (i in 1:M){
    cat("\n", rownames(object$BF.mat)[i], ":\n")
    cat("   call = \n")
    print(object$BF.call[[i]])
    cat("\n   log marginal likelihood = ", object$BF.logmarglike[i], "\n\n")
    
  }
  
}


"PostProbMod" <- function(BF, prior.probs=1){
  if (!is.BayesFactor(BF)){
    stop("BF is not of class BayesFactor\n")
  }

  M <- length(BF$BF.call)

  if (min(prior.probs) <= 0){
    stop("An element of prior.probs is non-positive\n") 
  }
  
  prior.probs <- rep(prior.probs, M)[1:M]
  prior.probs <- prior.probs / sum(prior.probs)
  
  lognumer <- BF$BF.logmarglike + log(prior.probs)
  maxlognumer <- max(lognumer)

  logpostprobs <- array(NA, M)
  denom <- 0
  for (i in 1:M){
    denom <- denom + exp(lognumer[i]-maxlognumer)
  }
  logdenom <- log(denom)
  
  for (i in 1:M){
    logpostprobs[i] <- (lognumer[i] - maxlognumer) - logdenom
  }
  postprobs <- exp(logpostprobs)

  names(postprobs) <- rownames(BF$BF.mat)

  return(postprobs)  
}
