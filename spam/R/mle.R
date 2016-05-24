# This is file ../spam/R/mle.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     



neg2loglikelihood.spam <- function(y, X, distmat, Covariance,
                                   beta, theta, Rstruct = NULL,...) {

  Sigma <- do.call(Covariance,list(distmat,theta))
  if (!is.spam(Sigma)){
    warning("'Covariance' should return a spam object. Forced to spam.")
    Sigma <- as.spam(Sigma)
  }
  if (is(Rstruct, "spam.chol.NgPeyton")) 
    cholS <- update.spam.chol.NgPeyton(Rstruct, Sigma, ...)
  else cholS <- chol.spam(Sigma, ...)

  n <- length(y)
  resid <- y-X%*%beta
  return( n * log(2*pi) +
         2*c(determinant.spam.chol.NgPeyton(cholS)$modulus) +
         sum(resid * solve.spam( cholS, resid))
         )
}


neg2loglikelihood <- function(y, X, distmat, Covariance,
                                   beta, theta, ...) {

  Sigma <- do.call(Covariance,list(distmat,theta))
  cholS <- chol(Sigma, ...)
  logdet <- sum(log(diag(cholS)))
  
  n <- length(y)
  resid <- y-X%*%beta
  
  return( n * log(2*pi) +
         2*logdet +
         sum(resid * backsolve(cholS, forwardsolve(cholS, resid,
            transpose=TRUE, upper.tri=TRUE),n))
         )
}





mle.spam <- function(y, X, distmat, Covariance,
                     beta0, theta0,
                     thetalower, thetaupper, optim.control=NULL,
                     Rstruct = NULL, hessian = FALSE,...) {
  
  
  if (!is(Rstruct, "spam.chol.NgPeyton")) {
    Sigma <- do.call(Covariance,
                     list(distmat,c(thetaupper[1],theta0[-1])))
    if (!is.spam(Sigma))
      stop("'Covariance' should return a spam object.")
   
    Rstruct <- chol.spam(Sigma, ...)
  }
  
  p <- dim(X)[2]
  n <- length(y)
    
  neg2loglikelihood <- function(fulltheta,...) {
    Sigma <- do.call(Covariance,list(distmat,fulltheta[-(1:p)]))
    cholS <- update.spam.chol.NgPeyton(Rstruct, Sigma, ...)

    resid <- y-X%*%fulltheta[1:p]
    return( n * log(2*pi) +
           2*c(determinant.spam.chol.NgPeyton(cholS)$modulus) +
           sum(resid * solve.spam( cholS, resid))
           )
  }

  return(optim(c(beta0,theta0),neg2loglikelihood,
               method = "L-BFGS-B",control = optim.control,
               lower=c(rep(-Inf,p),thetalower),
               upper=c(rep(Inf,p),thetaupper), hessian = hessian))
               
}

mle <- function(y, X, distmat, Covariance,
                beta0, theta0,
                thetalower, thetaupper, optim.control=NULL,
                hessian = FALSE,
                ...) {
    
  p <- dim(X)[2]
  n <- length(y)
    
  neg2loglikelihood <- function(fulltheta,...) {
    Sigma <- do.call(Covariance,list(distmat,fulltheta[-(1:p)]))
    cholS <- chol(Sigma, ...)
    logdet <- sum(log(diag(cholS)))

    resid <- y-X%*%fulltheta[1:p]
    return( n * log(2*pi) +
           2*logdet +
           sum(resid * backsolve(cholS, forwardsolve(cholS, resid,
              transpose=TRUE, upper.tri=TRUE),n))
           )
  }

  return(optim(c(beta0,theta0),neg2loglikelihood,
               method = "L-BFGS-B",control = optim.control,
               lower=c(rep(-Inf,p),thetalower),
               upper=c(rep(Inf,p),thetaupper), hessian = hessian))
               
}


mle.nomean.spam <- function(y, distmat, Covariance,
                     theta0,
                     thetalower, thetaupper, optim.control = NULL,
                     Rstruct = NULL, hessian = FALSE,...) {
  
  
  if (!is(Rstruct, "spam.chol.NgPeyton")) {
    Sigma <- do.call(Covariance,
                     list(distmat,c(thetaupper[1],theta0[-1])))
    if (!is.spam(Sigma))
      stop("'Covariance' should return a spam object.")
   
    Rstruct <- chol.spam(Sigma, ...)
  }
  
  n <- length(y)
    
  neg2loglikelihood <- function(theta,...) {
    Sigma <- do.call(Covariance,list(distmat,theta))
    cholS <- update.spam.chol.NgPeyton(Rstruct, Sigma, ...)

    return( n * log(2*pi) +
           2*c(determinant.spam.chol.NgPeyton(cholS)$modulus) +
           sum(y * solve.spam( cholS, y))
           )
  }

  return(optim(theta0,neg2loglikelihood,
               method = "L-BFGS-B",control = optim.control,
               lower=thetalower,    upper=thetaupper, hessian = hessian))
               
}




mle.nomean <- function(y, distmat, Covariance,
                       theta0,
                       thetalower, thetaupper, optim.control=NULL,
                       hessian = FALSE,
                       ...) {
    
  n <- length(y)
    
  neg2loglikelihood <- function(theta,...) {
    Sigma <- do.call(Covariance,list(distmat,theta))
    cholS <- chol(Sigma, ...)
    logdet <- sum(log(diag(cholS)))

    return( n * log(2*pi) +
           2*logdet +
           sum(y * backsolve(cholS, forwardsolve(cholS, y,
              transpose=TRUE, upper.tri=TRUE),n))
           )
  }

  return(optim(theta0,neg2loglikelihood,
               method = "L-BFGS-B",control = optim.control,
               lower=thetalower,    upper=thetaupper, hessian = hessian))
               
}
