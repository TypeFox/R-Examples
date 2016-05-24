##
##  PURPOSE:   (Log-)density of the multivariate Student t distribution
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  LOG:       20120124 created
##
##  FUNCTION:  dMVT
##
## ======================================================================

## *************************************************************
## dMVT
## *************************************************************
dMVT <- function(x, df, mu=0, Q=1, Sigma, log=FALSE)
{
  thispackage <- "mixAK"

  if (df <= 0) stop("df (degrees of freedom) must be positive")  
  
  ## dimension of the t distribution + precision matrix (lower triangle only)
  nx <- length(mu)
  if (nx == 1){    
    if (!missing(Sigma)) Q <- 1/Sigma
    if (length(Q) != 1) stop(paste("length(mu)=", nx, ", but length(Q)=", length(Q), "!", sep=""))    
  }else{
    if (missing(Sigma)){
      if (is.null(dim(Q))) stop(paste("length(mu)=", nx, ", but dim(Q) is NULL!", sep=""))
    }else{
      if (is.null(dim(Sigma))) stop(paste("length(mu)=", nx, ", but dim(Sigma) is NULL!", sep=""))
      Q <- chol2inv(chol(Sigma))
    }
    Q <- Q[lower.tri(Q, diag=TRUE)]
  }  
  
  ## number of points where to evaluate density
  if (is.null(dim(x))){
    if (nx == 1){
      npoints <- length(x)
    }
    else{
      if (length(x) != nx) stop(paste("Dimension of the MVT distribution is ", nx, " and x has length ", length(x), sep=""))
      npoints <- 1
    }  
  }
  else{
    if (ncol(x) != nx) stop(paste("Dimension of the MVT distribution is ", nx, " and x has ", ncol(x), " columns", sep=""))
    npoints <- nrow(x)
  }  

  DENSITY <- .C("dMVT1_R", value=double(npoints),
                           Q=as.double(Q),
                           work=double(nx),
                           err=integer(1),
                           x=as.double(t(x)),
                           unlog=as.integer(!log),
                           nu=as.double(df),
                           mu=as.double(mu),
                           nx=as.integer(nx),
                           npoints=as.integer(npoints),
                PACKAGE=thispackage)
  

  if (DENSITY$err) stop("Supplied precision/covariance matrix was not positive definite.")
  return(as.numeric(DENSITY$value))    
}  
