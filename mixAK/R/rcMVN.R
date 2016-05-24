##
##  PURPOSE:   Random number generation from the multivariate normal distribution
##             specified in a canonical form
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   05/11/2007
##
##  FUNCTION:  rcMVN
##
## ======================================================================

## *************************************************************
## rcMVN
## *************************************************************
rcMVN <- function(n, b=0, Q=1, Sigma)
{
  thispackage <- "mixAK"

  if (n <= 0) stop("n must be positive")
  
  ## dimension of the normal distribution + precision matrix (lower triangle only)
  nx <- length(b)
  if (nx == 1){    
    if (!missing(Sigma)) Q <- 1/Sigma
    if (length(Q) != 1) stop(paste("length(b)=", nx, ", but length(Q)=", length(Q), "!", sep=""))    
  }else{
    if (missing(Sigma)){
      if (is.null(dim(Q))) stop(paste("length(b)=", nx, ", but dim(Q) is NULL!", sep=""))
    }else{
      if (is.null(dim(Sigma))) stop(paste("length(b)=", nx, ", but dim(Q) is NULL!", sep=""))
      Q <- chol2inv(chol(Sigma))
    }
    Q <- Q[lower.tri(Q, diag=TRUE)]
  }  
  
  SAMPLE <- .C("rMVN2_R", x=double(nx*n),
                          mu=as.double(b),
                          log.dens=double(n),
                          Q=as.double(Q),
                          err=integer(1),
                          nx=as.integer(nx),
                          npoints=as.integer(n),
                PACKAGE=thispackage)

  if (SAMPLE$err) stop("Supplied precision/covariance matrix was not positive definite.")
  
  if (n == 1){
    names(SAMPLE$x) <- names(SAMPLE$mu) <- names(b)
    names(SAMPLE$log.dens) <- 1:n
  }else{
    SAMPLE$x <- matrix(SAMPLE$x, byrow=TRUE, ncol=nx, nrow=n)
    colnames(SAMPLE$x) <- names(SAMPLE$mu) <- names(b)
    names(SAMPLE$log.dens) <- rownames(SAMPLE$x) <- 1:n
  }
  
  return(list(x=SAMPLE$x, mean=SAMPLE$mu, log.dens=SAMPLE$log.dens))    
}  

