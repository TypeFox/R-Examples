##
##  PURPOSE:   Random number generation from the multivariate normal distribution
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   05/11/2007
##
##  FUNCTION:  rMVN
##
## ======================================================================

## *************************************************************
## rMVN
## *************************************************************
rMVN <- function(n, mean=0, Q=1, Sigma)
{
  thispackage <- "mixAK"

  if (n <= 0) stop("n must be positive")
  
  ## dimension of the normal distribution + precision matrix (lower triangle only)
  nx <- length(mean)
  if (nx == 1){    
    if (!missing(Sigma)) Q <- 1/Sigma
    if (length(Q) != 1) stop(paste("length(mean)=", nx, ", but length(Q)=", length(Q), "!", sep=""))    
  }else{
    if (missing(Sigma)){
      if (is.null(dim(Q))) stop(paste("length(mean)=", nx, ", but dim(Q) is NULL!", sep=""))
    }else{
      if (is.null(dim(Sigma))) stop(paste("length(mean)=", nx, ", but dim(Q) is NULL!", sep=""))
      Q <- chol2inv(chol(Sigma))
    }
    Q <- Q[lower.tri(Q, diag=TRUE)]
  }  

  if (any(mean != 0)) mu.nonZERO <- 1 else mu.nonZERO <- 0
  
  SAMPLE <- .C("rMVN1_R", x=double(nx*n),
                          log.dens=double(n),
                          Q=as.double(Q),
                          err=integer(1),
                          mu=as.double(mean),
                          nx=as.integer(nx),
                          mu.nonZERO=as.integer(mu.nonZERO),
                          npoints=as.integer(n),
                PACKAGE=thispackage)

  if (SAMPLE$err) stop("Supplied precision/covariance matrix was not positive definite.")
  
  if (n == 1){
    names(SAMPLE$x) <-  names(mean)
    names(SAMPLE$log.dens) <- 1:n
  }else{
    SAMPLE$x <- matrix(SAMPLE$x, byrow=TRUE, ncol=nx, nrow=n)
    colnames(SAMPLE$x) <- names(mean)
    names(SAMPLE$log.dens) <- rownames(SAMPLE$x) <- 1:n
  }
  
  return(list(x=SAMPLE$x, log.dens=SAMPLE$log.dens))    
}  



