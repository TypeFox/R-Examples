##
##  PURPOSE:   Computation of the pairwise joint densities
##             (plug-in version with supplied posterior summaries of mixture components)
##             * default method
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTION:  NMixPlugDensJoint2.default (28/05/2009) 
##
## ======================================================================

## *************************************************************
## NMixPlugDensJoint2.default
## *************************************************************
##
## Z ~ sum w[j] N(mu[j], Sigma[j])
## It computes bivariate joint densities of (X[d],X[c]), where
##    X[d] = scale$shift[d] + scale$scale[d] * Z[d]
##    X[c] = scale$shift[c] + scale$scale[c] * Z[c]
##
NMixPlugDensJoint2.default <- function(x, scale, w, mu, Sigma, ...)
{  
  ## Dimension of the normal mixture
  if (!is.list(x)) stop("x must be a list")
  p <- length(x)
  if (p < 2) stop("length of x must be 2 or more")
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

  ## Check mixture means and variances (REMEMBER, p >= 2)
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

  ## Adjust mixture means and variances with respect to scaling
  mu <- mu * matrix(rep(scale$scale, K), nrow=K, ncol=p, byrow=TRUE) + matrix(rep(scale$shift, K), nrow=K, ncol=p, byrow=TRUE)
  for (k in 1:K) Sigma[[k]] <- diag(scale$scale) %*% Sigma[[k]] %*% diag(scale$scale)

  ## Lengths of grids in each margin
  n <- sapply(x, length)
  if (any(n <= 0)) stop("incorrect x supplied")
  
  ## Compute pairwise joint densities
  RET <- list(x=x, dens=list())
  if (p == 2){
    GRID <- cbind(rep(x[[1]], n[2]), rep(x[[2]], each=n[1]))
    RET$dens[[1]] <- matrix(dMVNmixture2(x=GRID, weight=w, mean=mu, Sigma=Sigma), nrow=n[1], ncol=n[2])
    names(RET$dens) <- "1-2"
  }else{
    pp <- 1
    NAMEN <- character(0)
    for (m0 in 1:(p-1)){
      for (m1 in (m0+1):p){
        GRID <- cbind(rep(x[[m0]], n[m1]), rep(x[[m1]], each=n[m0]))
        MEAN <- mu[, c(m0, m1)]
        SIGMA <- list()
        for (k in 1:K) SIGMA[[k]] <- Sigma[[k]][c(m0, m1), c(m0, m1)]
        RET$dens[[pp]] <- matrix(dMVNmixture2(x=GRID, weight=w, mean=MEAN, Sigma=SIGMA), nrow=n[m0], ncol=n[m1])
        NAMEN <- c(NAMEN, paste(m0, "-", m1, sep=""))
        pp <- pp + 1
      }      
    }
    names(RET$dens) <- NAMEN    
  }  
    
  class(RET) <- "NMixPlugDensJoint2"
  return(RET)    
}  


