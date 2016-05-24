##
##  PURPOSE:   (Log-)density of the multivariate normal distribution
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   05/11/2007
##
##  FUNCTION:  dMVN
##
## ======================================================================

## *************************************************************
## dMVN
## *************************************************************
dMVN <- function(x, mean=0, Q=1, Sigma, log=FALSE)
{
  thispackage <- "mixAK"

  ## dimension of the normal distribution + precision matrix (lower triangle only)
  nx <- length(mean)
  if (nx == 1){    
    if (!missing(Sigma)) Q <- 1/Sigma
    if (length(Q) != 1) stop(paste("length(mean)=", nx, ", but length(Q)=", length(Q), "!", sep=""))    
  }else{
    if (missing(Sigma)){
      if (is.null(dim(Q))) stop(paste("length(mean)=", nx, ", but dim(Q) is NULL!", sep=""))
    }else{
      if (is.null(dim(Sigma))) stop(paste("length(mean)=", nx, ", but dim(Sigma) is NULL!", sep=""))
      Q <- chol2inv(chol(Sigma))
    }
    Q <- Q[lower.tri(Q, diag=TRUE)]
  }  

  if (any(mean != 0)) mu.nonZERO <- 1 else mu.nonZERO <- 0
  
  ## number of points where to evaluate density
  if (is.null(dim(x))){
    if (nx == 1){
      npoints <- length(x)
    }
    else{
      if (length(x) != nx) stop(paste("Dimension of the normal distribution is ", nx, " and x has length ", length(x), sep=""))
      npoints <- 1
    }  
  }
  else{
    if (ncol(x) != nx) stop(paste("Dimension of the normal distribution is ", nx, " and x has ", ncol(x), " columns", sep=""))
    npoints <- nrow(x)
  }  

  DENSITY <- .C("dMVN1_R", value=double(npoints),
                           Q=as.double(Q),
                           work=double(nx),
                           err=integer(1),
                           x=as.double(t(x)),
                           unlog=as.integer(!log),
                           mu=as.double(mean),
                           nx=as.integer(nx),
                           mu.nonZERO=as.integer(mu.nonZERO),
                           npoints=as.integer(npoints),
                PACKAGE=thispackage)
  

  if (DENSITY$err) stop("Supplied precision/covariance matrix was not positive definite.")
  return(as.numeric(DENSITY$value))    
}  
