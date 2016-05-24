##
##  PURPOSE:   Computation of the preditive univariate conditional cdf's (given one margin)
##             * default method
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   06/05/2010
##
##  FUNCTION:  NMixPredCondCDFMarg.default (06/05/2010) 
##
## ======================================================================

## *************************************************************
## NMixPredCondCDFMarg.default
## *************************************************************
##
## Z ~ sum w[j] N(mu[j], Sigma[j])
## It computes univariate conditional cdf's of X[d] | X[icond], where
##    X[d] = scale$shift[d] + scale$scale[d] * Z[d]
##

NMixPredCondCDFMarg.default <- function(x, icond, prob, scale, K, w, mu, Li, Krandom=FALSE, ...)
{
  thispackage <- "mixAK" 
  
  ## Dimension of the original normal mixture
  if (!is.list(x)) stop("x must be a list")
  p <- length(x)
  if (p <= 1) stop("length of x must be 2 or more")
  LTp <- p * (p + 1)/2

  if (p <= 1) stop("not sensible for univariate distribution")
  
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

  ## Check chains
  Kmax <- max(K)
  if (any(K <= 0)) stop("all K's must be positive")  
  if (Krandom){
    M <- length(K)
    sumK <- sum(K)
  }else{
    K <- K[1]
    M <- length(w)/K
    sumK <- K * M
  }
  if (length(w) != sumK) stop("incorrect w supplied")
  if (length(mu) != p*sumK) stop("incorrect mu supplied")  
  if (length(Li) != LTp*sumK) stop("incorrect Li supplied")  

  ## Lengths of grids in each margin
  n <- sapply(x, length)
  if (any(n <= 0)) stop("incorrect x supplied")

  ## Compute needed space
  ldens <- n[icond] + n[icond]*(sum(n[-icond]))

  ## Adjust grids with respect to scaling
  grid <- list()
  for (d in 1:p) grid[[d]] <- (x[[d]] - scale$shift[d])/scale$scale[d]

  ## Pointwise quantiles?
  if (missing(prob)){
    nquant <- 0
    prob <- 0
  }else{
    nquant <- length(prob)
    if (any(prob < 0)) stop("all prob values must be >= 0")
    if (any(prob > 1)) stop("all prob values must be <= 1")    
  }  
  
  if (Krandom) stop("not (yet) implemented for random K")

  ## Compute predictive densities
  RES <- .C("NMix_PredCondDensCDFMarg",
            dens      = double(ldens),
            qdens     = double(ifelse(nquant, ldens * nquant, 1)),
            err       = integer(1),
            calc_dens = as.integer(0),
            nquant    = as.integer(nquant),
            qprob     = as.double(prob),
            icond     = as.integer(icond-1),
            y         = as.double(unlist(grid)),
            p         = as.integer(p),
            n         = as.integer(n),
            chK       = as.integer(K),
            chw       = as.double(w),
            chmu      = as.double(mu),
            chLi      = as.double(Li),
            M         = as.integer(M),
            PACKAGE=thispackage)

  if (RES$err) stop("Something went wrong.")

  ## Lengths of grids other than x[[icond]]
  nmarg <- n[-icond]

  ## Number of values in RES$dens which preceed values for a specific margin
  voor.marg <- cumsum(c(n[icond], nmarg*n[icond]))
  voor.marg <- voor.marg[-length(voor.marg)]

  ## Add also 0 for conditioning margin and sort according indeces of margins
  margs <- (1:p)[-icond]
  voor <- rep(0, p)
  voor[margs] <- voor.marg

  ## Create resulting object
  RET <- list(x=x, icond=icond, cdf=list())
  for (t in 1:length(x[[icond]])){
    RET$cdf[[t]] <- list()

    for (m0 in (1:p)){
      if (m0 == icond){
        RET$cdf[[t]][[m0]] <- as.numeric(RES$dens[voor[m0] + t])/scale$scale[m0]
        next
      }  
      
      RET$cdf[[t]][[m0]] <- as.numeric(RES$dens[(voor[m0] + (t-1)*n[m0] + 1):(voor[m0] + t*n[m0])])
    }

    names(RET$cdf[[t]]) <- paste(1:p)    
  }

  ## Pointwise quantiles
  if (nquant){
    RET$prob <- prob
    for (i in 1:nquant){

      qnaam <- paste("q", prob[i]*100, "%", sep="")
      RET[[qnaam]] <- list()
      for (t in 1:length(x[[icond]])){
        RET[[qnaam]][[t]] <- list()

        for (m0 in (1:p)){
          if (m0 == icond){
            RET[[qnaam]][[t]][[m0]] <- as.numeric(RES$qdens[(i-1)*ldens + voor[m0] + t])/scale$scale[m0]
            next
          }  
      
          RET[[qnaam]][[t]][[m0]] <- as.numeric(RES$qdens[((i-1)*ldens + voor[m0] + (t-1)*n[m0] + 1):((i-1)*ldens + voor[m0] + t*n[m0])])
        }

        names(RET[[qnaam]][[t]]) <- paste(1:p)    
      }
      
    }      
  }  

  class(RET) <- "NMixPredCondCDFMarg"
  return(RET)    
}  


