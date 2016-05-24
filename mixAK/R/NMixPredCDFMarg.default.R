##
##  PURPOSE:   Computation of the predictive marginal (univariate) cumulative distribution functions
##             * default method
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/06/2009
##
##  FUNCTIONS: NMixPredCDFMarg.default (09/06/2009)
##             
## ======================================================================

## *************************************************************
## NMixPredCDFMarg.default
## *************************************************************
##
## Z ~ sum w[j] N(mu[j], (Li[j] %*% Li[j])^{-1})
## It computes marginal densities of X[d] = scale$shift[d] + scale$scale[d] * Z[d]
##
NMixPredCDFMarg.default <- function(x, scale, K, w, mu, Li, Krandom=TRUE, ...)
{  
  thispackage <- "mixAK"

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

  ## Total length of the grid
  lgrid <- sum(n)
  clgrids <- c(0, cumsum(n))

  ## Adjust grids with respect to scaling
  grid <- list()
  for (d in 1:p) grid[[d]] <- (x[[d]] - scale$shift[d])/scale$scale[d]
  
  RES <- .C("NMix_PredCDFMarg",
                cdf=double(lgrid),
                cdfK=double(lgrid*Kmax),
                freqK=integer(Kmax),
                propK=double(Kmax),
                dwork=double(LTp),
                err=integer(1),
                y=as.double(unlist(grid)),
                p=as.integer(p),
                n=as.integer(n),
                chK=as.integer(K),
                chw=as.double(w),
                chmu=as.double(mu),
                chLi=as.double(Li),
                M=as.integer(M),
                Kmax=as.integer(Kmax),
                Krandom=as.integer(Krandom),
            PACKAGE=thispackage)

  if (RES$err) stop("Something went wrong.")

  RET <- list(x=x, freqK=RES$freqK, propK=RES$propK, MCMC.length=M)
  if (p == 1){
    RET$cdf[[1]] <- RES$cdf
    RET$cdfK[[1]] <- list()
    for (j in 1:Kmax){  
      RET$cdfK[[1]][[j]] <- RES$cdfK[((j-1)*lgrid+1):(j*lgrid)]
    }
    names(RET$cdf) <- names(RET$cdfK) <- "1"
  }else{
    RET$cdf <- RET$cdfK <- list()
    for (m0 in 1:p){
      RET$cdf[[m0]] <- RES$cdf[(clgrids[m0]+1):clgrids[m0+1]]
      RET$cdfK[[m0]] <- list()
      for (j in 1:Kmax){
        RET$cdfK[[m0]][[j]] <- RES$cdfK[((j-1)*lgrid + clgrids[m0] + 1):((j-1)*lgrid + clgrids[m0+1])]
      }
    }
    names(RET$cdf) <- names(RET$cdfK) <- paste(1:p)
  }
  
  class(RET) <- "NMixPredCDFMarg"
  return(RET)    
}  


