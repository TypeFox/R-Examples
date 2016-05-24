##
##  PURPOSE:   Random number generation from the mixture of the multivariate normal distributions
##             * mixing performed in R
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   17/12/2007
##
##  FUNCTION:  rMVNmixture
##
## ======================================================================

## *************************************************************
## rMVNmixture
## *************************************************************
rMVNmixture <- function(n, weight, mean, Q, Sigma)
{
  thispackage <- "mixAK"

  if (n <= 0) stop("n must be positive")
  
  ## number of components of the mixture + checking the weights
  if (any(weight < 0)) stop("weights must be non-negative")  
  K <- length(weight)
  if (K == 0) stop("weight is of zero length")
  weight <- weight/sum(weight)

  ## dimension of the normal distribution + precision/covariance matrix
  UNIVARIATE <- FALSE
  if (!is.matrix(mean)) UNIVARIATE <- TRUE
  else                  if (ncol(mean) == 1) UNIVARIATE <- TRUE    
  if (UNIVARIATE){      ### univariate
    p <- 1
    if (length(mean) != K) stop(paste("mean must be of length ", K, sep=""))
    
    if (missing(Sigma)){
      if (missing(Q)) stop("Sigma or Q must be given")
      if (length(Q) != K) stop(paste("Q must be of length ", K, sep=""))      
      if (is.list(Q)){
        lQ <- sapply(Q, length)
        if (any(lQ != 1)) stop("all Q elements must be of length 1")
        Q <- unlist(Q)        
      }
      if (any(Q <= 0)) stop("all Q elements must be positive")
      degener <- is.infinite(Q)
      Sigma <- 1/Q
    }else{
      if (length(Sigma) != K) stop(paste("Sigma must be of length ", K, sep=""))      
      if (is.list(Sigma)){
        lSigma <- sapply(Sigma, length)
        if (any(lSigma != 1)) stop("all Sigma elements must be of length 1")
        Sigma <- unlist(Sigma)        
      }
      if (any(Sigma < 0)) stop("all Sigma elements must be nonnegative")
      degener <- (Sigma==0)
      Q <- 1/Sigma      
    }  
  }else{                       ### multivariate
    p <- ncol(mean)
    if (nrow(mean) != K) stop(paste("mean must have ", K, " rows", sep=""))

    if (missing(Sigma)){
      if (missing(Q)) stop("Sigma or Q must be given")            
      if (is.matrix(Q)){
        if (K != 1) stop("Q must be a list of matrices")
        Q <- list(Q)
      }
      if (length(Q) != K) stop(paste("Q must be of length ", K, sep=""))
      Sigma <- QLT <- list()      
      for (j in 1:K){
        if (!is.matrix(Q[[j]])) stop("all elements of Q must be matrices")
        if (nrow(Q[[j]]) != p | ncol(Q[[j]]) != p) stop(paste("all elements of Q must be squared matrices with ", p, " rows and columns", sep=""))
        Sigma[[j]] <- chol2inv(chol(Q[[j]]))
        QLT[[j]] <- Q[[j]][lower.tri(Q[[j]], diag=TRUE)]
      }                  
    }else{
      if (is.matrix(Sigma)){
        if (K != 1) stop("Sigma must be a list of matrices")
        Sigma <- list(Sigma)
      }
      if (length(Sigma) != K) stop(paste("Sigma must be of length ", K, sep=""))
      Q <- QLT <- list()      
      for (j in 1:K){
        if (!is.matrix(Sigma[[j]])) stop("all elements of Sigma must be matrices")
        if (nrow(Sigma[[j]]) != p | ncol(Sigma[[j]]) != p) stop(paste("all elements of Sigma must be squared matrices with ", p, " rows and columns", sep=""))
        Q[[j]] <- chol2inv(chol(Sigma[[j]]))
        QLT[[j]] <- Q[[j]][lower.tri(Q[[j]], diag=TRUE)]
      }                  
    }     
  }  

  ## sample the components
  r <- sample(1:K, size=n, replace=TRUE, prob=weight)
  
  ## sample the points from the components
  if (p == 1){
    StdDev <- sqrt(Sigma)
    x <- rnorm(n, mean=mean[r], sd=StdDev[r])    
  }else{
    x <- matrix(NA, nrow=n, ncol=p)
    for (j in 1:K){
      Nj <- sum(r == j)
      if (Nj > 0){
        x[r==j,] <- matrix(
                       .C("rMVN1_R", x=double(p*Nj),
                                     log.dens=double(Nj),
                                     Q=as.double(QLT[[j]]),
                                     work=double(p),
                                     err=integer(1),
                                     mu=as.double(mean[j,]),
                                     nx=as.integer(p),
                                     mu.nonZERO=as.integer(any(mean[j,] != 0)),
                                     npoints=as.integer(Nj),
                          PACKAGE=thispackage)$x,
                          nrow=Nj, ncol=p, byrow=TRUE)
      }  
    }
    if (n == 1) x <- as.numeric(x)
  }  
  
  return(x)    
}
