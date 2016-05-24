##
##  PURPOSE:   (Log-)density of the mixture of the multivariate normal distributions
##             * mixing performed in R
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   17/12/2007
##
##  FUNCTION:  dMVNmixture
##
## ======================================================================

## *************************************************************
## dMVNmixture
## *************************************************************
dMVNmixture <- function(x, weight, mean, Q, Sigma, log=FALSE)
{
  thispackage <- "mixAK"

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

  ## number of points where to evaluate density
  if (!is.matrix(x)){
    if (p == 1){
      npoints <- length(x)
    }else{
      if (length(x) != p) stop(paste("Dimension of the normal distribution is ", p, " and x has length ", length(x), sep=""))
      npoints <- 1
    }  
  }else{
    if (ncol(x) != p) stop(paste("Dimension of the normal distribution is ", p, " and x has ", ncol(x), " columns", sep=""))
    npoints <- nrow(x)
  }  

  ## evaluate the density
  DENSITY <- matrix(0, nrow=npoints, ncol=K)  
  if (p == 1){
    for (j in 1:K){
      if (degener[j]) DENSITY[,j] <- weight[j]*(x==mean[j])
      else            DENSITY[,j] <- weight[j]*dnorm(x, mean=mean[j], sd=sqrt(Sigma[j]))      
    }
  }else{    
    for (j in 1:K){      
      DENSITY[,j] <- weight[j]*.C("dMVN1_R", value=double(npoints),
                                             Q=as.double(QLT[[j]]),
                                             work=double(p),
                                             err=integer(1),
                                             x=as.double(t(x)),
                                             unlog=as.integer(1),
                                             mu=as.double(mean[j,]),
                                             nx=as.integer(p),
                                             mu.nonZERO=as.integer(any(mean[j,] != 0)),
                                             npoints=as.integer(npoints),
                                  PACKAGE=thispackage)$value      
    }
  }  

  DENSITY <- apply(DENSITY, 1, sum)
  if (log) DENSITY <- log(DENSITY)
  
  return(DENSITY)  
}

