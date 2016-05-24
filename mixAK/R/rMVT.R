##
##  PURPOSE:   Random number generation from the multivariate Student t distribution
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  LOG:       20111208  created
##
##  FUNCTION:  rMVT
##
## ======================================================================

## *************************************************************
## rMVT
## *************************************************************
rMVT <- function(n, df, mu=0, Q=1, Sigma)
{
  thispackage <- "mixAK"

  if (n <= 0) stop("n must be positive")
  if (df <= 0) stop("df (degrees of freedom) must be positive")  
  
  ## dimension of the MVT distribution + precision matrix (lower triangle only)
  nx <- length(mu)
  if (nx == 1){    
    if (!missing(Sigma)) Q <- 1/Sigma
    if (length(Q) != 1) stop(paste("length(mu)=", nx, ", but length(Q)=", length(Q), "!", sep=""))    
  }else{
    if (missing(Sigma)){
      if (is.null(dim(Q))) stop(paste("length(mu)=", nx, ", but dim(Q) is NULL!", sep=""))
    }else{
      if (is.null(dim(Sigma))) stop(paste("length(mu)=", nx, ", but dim(Q) is NULL!", sep=""))
      Q <- chol2inv(chol(Sigma))
    }
    Q <- Q[lower.tri(Q, diag=TRUE)]
  }  
  
  SAMPLE <- .C("rMVT1_R", x=double(nx*n),
                          log.dens=double(n),
                          Q=as.double(Q),
                          err=integer(1),
                          nu=as.double(df),
                          mu=as.double(mu),
                          nx=as.integer(nx),
                          npoints=as.integer(n),
                PACKAGE=thispackage)

  if (SAMPLE$err) stop("Supplied precision/covariance matrix was not positive definite.")
  
  if (n == 1){
    names(SAMPLE$x) <-  names(mu)
    names(SAMPLE$log.dens) <- 1:n
  }else{
    SAMPLE$x <- matrix(SAMPLE$x, byrow=TRUE, ncol=nx, nrow=n)
    colnames(SAMPLE$x) <- names(mu)
    names(SAMPLE$log.dens) <- rownames(SAMPLE$x) <- 1:n
  }
  
  return(list(x=SAMPLE$x, log.dens=SAMPLE$log.dens))    
}  
