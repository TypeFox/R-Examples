##
##  PURPOSE:   Computation of the marginal densities
##             (plug-in version with supplied posterior summaries of mixture components)
##             * default method
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTION:  NMixPlugDensMarg.default (28/05/2009) 
##
## ======================================================================

## *************************************************************
## NMixPlugDensMarg.default
## *************************************************************
##
## Z ~ sum w[j] N(mu[j], Sigma[j])
## It computes marginal densities of X[d], where
##    X[d] = scale$shift[d] + scale$scale[d] * Z[d]
##
NMixPlugDensMarg.default <- function(x, scale, w, mu, Sigma, ...)
{  
  ## Dimension of the normal mixture
  if (is.numeric(x)) x <- list(x1=x)  
  if (!is.list(x)) stop("x must be a list")
  p <- length(x)
  if (p < 1) stop("length of x must be 1 or more")
  LTp <- p * (p + 1)/2

  if (is.null(names(x))) names(x) <- paste("x", (1:p), sep="")
  
  ## scale
  if (missing(scale)) scale <- list(shift=rep(0, p), scale=rep(1, p))
  if (!is.list(scale)) stop("scale must be a list")
  if (length(scale) != 2) stop("scale must have 2 components")
  inscale <- names(scale)  
  iscale.shift <- match("shift", inscale, nomatch=NA)
  iscale.scale <- match("scale", inscale, nomatch=NA)
  if (is.na(iscale.shift)) stop("scale$shift is missing")
  if (length(scale$shift) == 1) scale$shift <- rep(scale$shift, p)
  if (length(scale$shift) != p) stop(paste("scale$shift must be a vector of length ", p, sep=""))    
  if (is.na(iscale.scale)) stop("scale$scale is missing")
  if (length(scale$scale) == 1) scale$scale <- rep(scale$scale, p)
  if (length(scale$scale) != p) stop(paste("scale$scale must be a vector of length ", p, sep=""))
  if (any(scale$scale <= 0)) stop("all elements of scale$scale must be positive")

  ## number of mixture components
  K <- length(w)

  ## Check mixture weights
  if (any(w < 0) | any(w > 1)) stop("weights must lie between zero and 1")
  if (abs(sum(w) - 1) > 1e-5) warning("sum of weights differs from 1")

  ## Check mixture means and variances
  ## (make them numeric if p=1)
  ## Adjust means and variances for scaling
  if (p == 1){

    ### Check
    if (length(mu) != K) stop("incorrect mu")
    mu <- as.numeric(mu)
    if (is.list(Sigma)){
      if (any(sapply(Sigma, length) != 1)) stop("incorrect Sigma")
      Sigma <- unlist(Sigma)
    }
    if (length(Sigma) != K) stop("incorrect Sigma")

    ### Scale adjustment
    mu <- mu * scale$scale + scale$shift
    Sigma <- Sigma * scale$scale^2
    
  }else{

    ### Check
    if (K == 1){
      if (length(mu) != p) stop("incorrect mu")
      mu <- matrix(mu, nrow=K, ncol=p)
      if (!is.list(Sigma)) Sigma <- list(Sigma)    
    }
    if (nrow(mu) != K) stop(paste("mu must have ", K, " rows", sep=""))
    if (ncol(mu) != p) stop(paste("mu must have ", p, " columns", sep=""))
    if (any(!sapply(Sigma, is.matrix))) stop("all Sigma's must be matrices")
    if (any(sapply(Sigma, nrow) != p)) stop(paste("all Sigma's must have ", p, " rows", sep=""))
    if (any(sapply(Sigma, ncol) != p)) stop(paste("all Sigma's must have ", p, " columns", sep=""))

    ### Scale adjustment
    mu <- mu * matrix(rep(scale$scale, K), nrow=K, ncol=p, byrow=TRUE) + matrix(rep(scale$shift, K), nrow=K, ncol=p, byrow=TRUE)
    for (k in 1:K) Sigma[[k]] <- diag(scale$scale) %*% Sigma[[k]] %*% diag(scale$scale)    
  }
      
  ## Lengths of grids in each margin
  n <- sapply(x, length)
  if (any(n <= 0)) stop("incorrect x supplied")
  
  ## Compute marginal densities
  RET <- list(x=x, dens=list())
  if (p == 1){
    RET$dens[[1]] <- dMVNmixture2(x=x[[1]], weight=w, mean=mu, Sigma=Sigma)
  }else{
    for (m0 in 1:p){
      MEAN <- as.numeric(mu[,m0])
      SIGMA <- Sigma[[1]][m0, m0]
      if (K >= 2) for (k in 2:K) SIGMA <- c(SIGMA, Sigma[[k]][m0, m0])
      RET$dens[[m0]] <- dMVNmixture2(x=x[[m0]], weight=w, mean=MEAN, Sigma=SIGMA)
    }  
  }    
  names(RET$dens) <- paste(1:p)
    
  class(RET) <- "NMixPlugDensMarg"
  return(RET)    
}  


