##
##  PURPOSE:   Computation of the pairwise bivariate conditional densities (given one margin)
##             (plug-in version with supplied posterior summaries of mixture components)
##             * default method
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   29/05/2009
##
##  FUNCTION:  NMixPlugCondDensJoint2.default (29/05/2009) 
##
## ======================================================================

## *************************************************************
## NMixPlugCondDensJoint2.default
## *************************************************************
##
## Z ~ sum w[j] N(mu[j], Sigma[j])
## It computes pairwise conditional densities of (X[c], X[d]) | X[icond], where
##    X[i] = scale$shift[i] + scale$scale[i] * Z[i]
##

NMixPlugCondDensJoint2.default <- function(x, icond, scale, w, mu, Sigma, ...)
{
  ## Dimension of the normal mixture
  if (!is.list(x)) stop("x must be a list")
  p <- length(x)
  if (p <= 2) stop("length of x must be 3 or more")
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

  ## Compute conditional densities for remaining pairs of margins
  RET <- list(x=x, icond=icond, dens=list())
  for (t in 1:length(x[[icond]])){
    RET$dens[[t]] <- list()
  
    pp <- 1
    NAMEN <- character(0)
    for (m0 in 1:(p-1)){      
      for (m1 in (m0+1):p){
        if (m0 == icond | m1 ==icond){
          RET$dens[[t]][[pp]] <- NA
          NAMEN <- c(NAMEN, paste(m0, "-", m1, sep=""))          
          pp <- pp + 1
          next
        }  
        
        MEAN <- mu[, c(m0, m1, icond)]
        SIGMA <- list()
        for (k in 1:K) SIGMA[[k]] <- Sigma[[k]][c(m0, m1, icond), c(m0, m1, icond)]
        GRID <- cbind(rep(x[[m0]], n[m1]), rep(x[[m1]], each=n[m0]), rep(x[[icond]][t], n[m0]*n[m1]))
        lognumer <- dMVNmixture2(x=GRID, weight=w, mean=MEAN, Sigma=SIGMA, log=TRUE)        
        RET$dens[[t]][[pp]] <- matrix(exp(lognumer - logdenom[t]), nrow=n[m0], ncol=n[m1])
        NAMEN <- c(NAMEN, paste(m0, "-", m1, sep=""))
        pp <- pp + 1
      }      
    }
    names(RET$dens[[t]]) <- NAMEN
  }  
    
  class(RET) <- "NMixPlugCondDensJoint2"
  return(RET)    
}
