##
##  PURPOSE:   Computation of the predictive pairwise joint densities
##             * default method
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   03/12/2007
##
##  FUNCTION:  NMixPredDensJoint2.default (03/12/2007) 
##
## ======================================================================

## *************************************************************
## NMixPredDensJoint2.default
## *************************************************************
##
## Z ~ sum w[j] N(mu[j], (Li[j] %*% Li[j])^{-1})
## It computes bivariate joint densities of (X[d],X[c]), where
##    X[d] = scale$shift[d] + scale$scale[d] * Z[d]
##    X[c] = scale$shift[c] + scale$scale[c] * Z[c]
##
NMixPredDensJoint2.default <- function(x, scale, K, w, mu, Li, Krandom=TRUE, ...)
{  
  thispackage <- "mixAK"

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
    for (m1 in (m0+1):p){
      lgrids <- c(lgrids, n[m0]*n[m1])
    }
  }

  ## Total length of the grid
  lgrid <- sum(lgrids)
  clgrids <- c(0, cumsum(lgrids))

  ## Adjust grids with respect to scaling and compute det(S), det(S) is not needed here
  grid <- list()
  for (d in 1:p) grid[[d]] <- (x[[d]] - scale$shift[d])/scale$scale[d]
  #detS <- prod(scale$scale)
  
  RES <- .C("NMix_PredDensJoint2",
                dens=double(lgrid),
                densK=double(lgrid*Kmax),
                freqK=integer(Kmax),
                propK=double(Kmax),
                dwork=double(6+LTp+2+3),
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
  if (p == 2){
    RET$dens[[1]] <- matrix(RES$dens, nrow=n[1], ncol=n[2]) / (scale$scale[1]*scale$scale[2])
    RET$densK[[1]] <- list()
    for (j in 1:Kmax){  
      RET$densK[[1]][[j]] <- matrix(RES$densK[((j-1)*lgrid+1):(j*lgrid)], nrow=n[1], ncol=n[2]) / (scale$scale[1]*scale$scale[2])
    }
    names(RET$dens) <- names(RET$densK) <- "1-2"
  }else{
    RET$dens <- RET$densK <- list()
    pp <- 1
    NAMEN <- character(0)
    for (m0 in 1:(p-1)){
      for (m1 in (m0+1):p){
        RET$dens[[pp]] <- matrix(RES$dens[(clgrids[pp]+1):clgrids[pp+1]], nrow=n[m0], ncol=n[m1]) / (scale$scale[m0]*scale$scale[m1])
        RET$densK[[pp]] <- list()
        NAMEN <- c(NAMEN, paste(m0, "-", m1, sep=""))
        for (j in 1:Kmax){
          RET$densK[[pp]][[j]] <- matrix(RES$densK[((j-1)*lgrid + clgrids[pp] + 1):((j-1)*lgrid + clgrids[pp+1])], nrow=n[m0], ncol=n[m1]) / (scale$scale[m0]*scale$scale[m1])
        }
        pp <- pp + 1
      }      
    }
    names(RET$dens) <- names(RET$densK) <- NAMEN    
  }
  
  class(RET) <- "NMixPredDensJoint2"
  return(RET)    
}  


