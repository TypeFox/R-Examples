##
##  PURPOSE:   Computation of the preditive pairwise bivariate conditional densities (given one margin)
##             * default method
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   01/06/2009
##
##  FUNCTION:  NMixPredCondDensJoint2.default (01/06/2009) 
##
## ======================================================================

## *************************************************************
## NMixPredCondDensJoint2.default
## *************************************************************
##
## Z ~ sum w[j] N(mu[j], Sigma[j])
## It computes univariate conditional densities of (X[d], X[c]) | X[icond], where
##    X[d] = scale$shift[d] + scale$scale[d] * Z[d]
##

NMixPredCondDensJoint2.default <- function(x, icond, scale, K, w, mu, Li, Krandom=FALSE, ...)
{
  thispackage <- "mixAK" 
  
  ## Dimension of the original normal mixture
  if (!is.list(x)) stop("x must be a list")
  p <- length(x)
  if (p <= 1) stop("length of x must be 2 or more")
  LTp <- p * (p + 1)/2

  if (p <= 2) stop("not sensible for uni- or bivariate distribution")
  
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

  ## Lengths of grids for each pair
  lgrids <- numeric(0)
  for (m0 in 1:(p-1)){
    if (m0 == icond) next
    for (m1 in (m0+1):p){
      if (m1 == icond) next
      lgrids <- c(lgrids, n[m0]*n[m1])
    }
  }  
  
  ## Compute needed space
  ldens <- n[icond] + n[icond]*(sum(lgrids))
  lwork <- 3 + LTp + ldens

  ## Adjust grids with respect to scaling
  grid <- list()
  for (d in 1:p) grid[[d]] <- (x[[d]] - scale$shift[d])/scale$scale[d]

  if (Krandom) stop("not (yet) implemented for random K")
  
  ## Compute predictive densities
  RES <- .C("NMix_PredCondDensJoint2",
                dens=double(ldens),
                dwork=double(lwork),
                err=integer(1),
                icond=as.integer(icond-1),
                y=as.double(unlist(grid)),
                p=as.integer(p),
                n=as.integer(n),
                chK=as.integer(K),
                chw=as.double(w),
                chmu=as.double(mu),
                chLi=as.double(Li),
                M=as.integer(M),
            PACKAGE=thispackage)

  if (RES$err) stop("Something went wrong.")
  
  ## Number of values in RES$dens which preceed values for a specific pair of margins
  voor.pair <- c(cumsum(c(n[icond], lgrids*n[icond])))
  voor.pair <- voor.pair[-length(voor.pair)]
  
  ## Create resulting object
  RET <- list(x=x, icond=icond, dens=list())
  for (t in 1:length(x[[icond]])){
    RET$dens[[t]] <- list()

    pp <- 1
    ppair <- 1
    NAMEN <- character(0)
    for (m0 in (1:(p-1))){
      for (m1 in ((m0+1):p)){
        if (m0 == icond | m1 ==icond){
          RET$dens[[t]][[pp]] <- NA
          NAMEN <- c(NAMEN, paste(m0, "-", m1, sep=""))          
          pp <- pp + 1
          next
        }

        RET$dens[[t]][[pp]] <- matrix(RES$dens[(voor.pair[ppair] + (t-1)*lgrids[ppair] + 1):(voor.pair[ppair] + t*lgrids[ppair])], nrow=n[m0], ncol=n[m1], byrow=TRUE) / (scale$scale[m0]*scale$scale[m1])
        NAMEN <- c(NAMEN, paste(m0, "-", m1, sep=""))
        pp <- pp + 1
        ppair <- ppair + 1
      }  
    }  
    names(RET$dens[[t]]) <- NAMEN    
  }  

  class(RET) <- "NMixPredCondDensJoint2"
  return(RET)    
}  


