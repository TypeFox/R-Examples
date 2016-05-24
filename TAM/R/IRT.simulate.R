
#####################################################
# S3 method
IRT.simulate <- function (object, ...) {
  UseMethod("IRT.simulate")
}
#####################################################


######################################################
simulate_mml <- function(object, iIndex = NULL, theta = NULL, nobs = NULL, ...){
  
  nnodes <- nobs
  A <- object$A
  xsi <- object$xsi$xsi
  B <- object$B
  ndim <- dim(B)[3]
  guess <- object$guess
  maxK <- dim(A)[2]
  
  if(is.null(iIndex)){
    iIndex <- 1:dim(A)[1]
  }
  nI <- length(iIndex)
  
  if(is.null(theta)){  
    if(is.null(nnodes)){
      nnodes <- nrow(object$person)  
    }
    
    t.mean <- object$beta
    t.sigma <- object$variance
    if(ndim == 1){
      theta <- stats::rnorm(nnodes, mean = t.mean, sd = sqrt(t.sigma))  
    } else {
      theta <- MASS::mvrnorm(nnodes, mu = t.mean, Sigma = t.sigma)
    }  
    theta <-  matrix(theta, ncol = ndim)
  }
  
  if(is.null(dim(theta)) & ndim == 1){
    warning("Theta points should be either NULL or a matrix of same dimensionality as the trait distribution.")
    theta <- matrix(theta, ncol = ndim)
    
  }  
  if(ncol(theta) != ndim){
    warning("Theta points should be either NULL or a matrix of same dimensionality as the trait distribution.")
    theta <- matrix(theta, ncol = ndim)
  }
  nnodes <- nrow(theta)
  
  #****
  # calculate probs
  p <- IRT.irfprob(object = class(object), A = A, B = B, xsi = xsi, theta = theta,
                   guess = guess, nnodes = nnodes, iIndex = iIndex, maxK = maxK, ... )
  
  #****
  # simulate data
  res <- matrix( stats::runif(nnodes * nI), nrow = nnodes, ncol = nI)
  for(ii in 1:nI){
    cat.success.ii <- (res[, ii] > t(apply(p[ii, , ], 2, cumsum)))
    res[, ii] <-  c(cat.success.ii %*% rep(1, maxK))
  }
  
  res[ is.na(object$resp) ] <- NA
  
  #****
  # Output
  class(res) <- "IRT.simulate"
  return(res)
}

######################################################################
# S3 methods
IRT.simulate.tam.mml <- simulate_mml
IRT.simulate.tam.mml.2pl <- simulate_mml
IRT.simulate.tam.mml.3pl <- simulate_mml
IRT.simulate.tam.mml.mfr <- simulate_mml
######################################################################