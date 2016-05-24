##
##  PURPOSE:   Random number generation from the mixture of the multivariate normal distributions
##             * mixing performed in C++ code
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   07/11/2008
##
##  FUNCTION:  rMVNmixture2
##
## ======================================================================

## *************************************************************
## rMVNmixture2
## *************************************************************
rMVNmixture2 <- function(n, weight, mean, Q, Sigma)
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
      Sigma <- list()
      QLT <- numeric(0)
      for (j in 1:K){
        if (!is.matrix(Q[[j]])) stop("all elements of Q must be matrices")
        if (nrow(Q[[j]]) != p | ncol(Q[[j]]) != p) stop(paste("all elements of Q must be squared matrices with ", p, " rows and columns", sep=""))
        Sigma[[j]] <- chol2inv(chol(Q[[j]]))
        QLT <- c(QLT, Q[[j]][lower.tri(Q[[j]], diag=TRUE)])
      }                  
    }else{
      if (is.matrix(Sigma)){
        if (K != 1) stop("Sigma must be a list of matrices")
        Sigma <- list(Sigma)
      }
      if (length(Sigma) != K) stop(paste("Sigma must be of length ", K, sep=""))
      Q <- list()
      QLT <- numeric(0)
      for (j in 1:K){
        if (!is.matrix(Sigma[[j]])) stop("all elements of Sigma must be matrices")
        if (nrow(Sigma[[j]]) != p | ncol(Sigma[[j]]) != p) stop(paste("all elements of Sigma must be squared matrices with ", p, " rows and columns", sep=""))
        Q[[j]] <- chol2inv(chol(Sigma[[j]]))
        QLT <- c(QLT, Q[[j]][lower.tri(Q[[j]], diag=TRUE)])
      }                  
    }     
  }  

  ## sample
  if (p == 1){
    SAMPLE <- .C("rmixNorm_R", x      = double(n),
                               dens   = double(n),
                               cumw   = double(K),
                               K      = as.integer(K),
                               w      = as.double(weight),
                               mu     = as.double(mean),
                               sigma  = as.double(sqrt(Sigma)),
                               npoints= as.integer(n),
                  PACKAGE=thispackage)
    x <- SAMPLE$x    
  }else{
    SAMPLE <- .C("rmixMVN_R", x      = double(p*n),
                              dens   = double(n),
                              w.dets = as.double(weight),
                              cumw   = double(K),
                              Li     = as.double(QLT),
                              work   = double(p),
                              err    = integer(1),
                              K      = as.integer(K),
                              mu     = as.double(t(mean)),
                              nx     = as.integer(p),
                              npoints= as.integer(n),
                  PACKAGE=thispackage)
    if (SAMPLE$err) stop("Something went wrong.")
    x <- matrix(SAMPLE$x, nrow=n, ncol=p, byrow=TRUE)
  }

  return(list(x=x, dens=SAMPLE$dens))    
}
