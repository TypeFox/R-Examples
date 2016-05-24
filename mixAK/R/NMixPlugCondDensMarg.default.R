##
##  PURPOSE:   Computation of the univariate conditional densities (given one margin)
##             (plug-in version with supplied posterior summaries of mixture components)
##             * default method
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTION:  NMixPlugCondDensMarg.default (28/05/2009) 
##
## ======================================================================

## *************************************************************
## NMixPlugCondDensMarg.default
## *************************************************************
##
## Z ~ sum w[j] N(mu[j], Sigma[j])
## It computes univariate conditional densities of X[d] | X[icond], where
##    X[d] = scale$shift[d] + scale$scale[d] * Z[d]
##

NMixPlugCondDensMarg.default <- function(x, icond, scale, w, mu, Sigma, ...)
{  
  ## Dimension of the original normal mixture
  if (!is.list(x)) stop("x must be a list")
  p <- length(x)
  if (p <= 1) stop("length of x must be 2 or more")
  LTp <- p * (p + 1)/2

  if (icond < 1 | icond > p) stop(paste("icond must be between 1 and ", p, sep=""))  
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
  if (nrow(mu) != K) stop(paste("mu must have ", K, " rows", sep=""))
  if (ncol(mu) != p) stop(paste("mu must have ", p, " columns", sep=""))
  if (any(!sapply(Sigma, is.matrix))) stop("all Sigma's must be matrices")
  if (any(sapply(Sigma, nrow) != p)) stop(paste("all Sigma's must have ", p, " rows", sep=""))
  if (any(sapply(Sigma, ncol) != p)) stop(paste("all Sigma's must have ", p, " columns", sep=""))

  ## Adjust means and variances for scaling
  mu <- mu * matrix(rep(scale$scale, K), nrow=K, ncol=p, byrow=TRUE) + matrix(rep(scale$shift, K), nrow=K, ncol=p, byrow=TRUE)
  for (k in 1:K) Sigma[[k]] <- diag(scale$scale) %*% Sigma[[k]] %*% diag(scale$scale)    
     
  ## Lengths of grids in each margin
  n <- sapply(x, length)
  if (any(n <= 0)) stop("incorrect x supplied")

  ## Compute marginal log-density for denominator of the conditional densities
  MEAN <- as.numeric(mu[,icond])
  SIGMA <- Sigma[[1]][icond, icond]
  if (K >= 2) for (k in 2:K) SIGMA <- c(SIGMA, Sigma[[k]][icond, icond])
  logdenom <- dMVNmixture2(x=x[[icond]], weight=w, mean=MEAN, Sigma=SIGMA, log=TRUE)    
  
  ## Compute conditional densities for remaining margins
  RET <- list(x=x, icond=icond, dens=list())
  for (t in 1:length(x[[icond]])){
    RET$dens[[t]] <- list()

    for (m0 in (1:p)){

      if (m0 == icond){
        RET$dens[[t]][[m0]] <- exp(logdenom[t])
        next
      }  
      
      MEAN <- mu[, c(m0, icond)]
      SIGMA <- list()
      for (k in 1:K) SIGMA[[k]] <- Sigma[[k]][c(m0, icond), c(m0, icond)]
      GRID <- cbind(x[[m0]], rep(x[[icond]][t], length(x[[m0]])))
      lognumer <- dMVNmixture2(x=GRID, weight=w, mean=MEAN, Sigma=SIGMA, log=TRUE)
      RET$dens[[t]][[m0]] <- exp(lognumer - logdenom[t])
    }

    names(RET$dens[[t]]) <- paste(1:p)    
  }  
    
  class(RET) <- "NMixPlugCondDensMarg"
  return(RET)    
}  


