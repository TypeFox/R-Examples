##
##  PURPOSE:   Computation of the predictive conditional pairwise densities given one margin
##              * method for objects of class NMixMCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   01/06/2009
##
##  FUNCTIONS: NMixPredCondDensJoint2.NMixMCMC (01/06/2009)
##             
## ======================================================================

## *************************************************************
## NMixPredCondDensJoint2.NMixMCMC
## *************************************************************
NMixPredCondDensJoint2.NMixMCMC <- function(x, icond, grid, lgrid=50, scaled=FALSE, ...)
{
  if (missing(icond)) stop("icond must be given")
  if (x$dim <= 2) stop("not applicable for uni- or bivariate mixtures")  

  if (x$nx_w > 1) stop("This function has not (yet) been implemented if a factor covariate on mixture weights is present.")
  
  if (missing(grid)){
    grid <- list()
    if (scaled){
      for (i in 1:x$dim){
        rangeGrid <- x$summ.z.Mean["Median", i] + c(-3.5, 3.5)*x$summ.z.SDCorr["Median", (i-1)*(2*x$dim - i + 2)/2 + 1]
        grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
      }
    }else{
      for (i in 1:x$dim){
        rangeGrid <- x$summ.y.Mean["Median", i] + c(-3.5, 3.5)*x$summ.y.SDCorr["Median", (i-1)*(2*x$dim - i + 2)/2 + 1]
        grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
      }
    }
    names(grid) <- paste("x", 1:x$dim, sep="")    
  }

  if (!is.list(grid)) stop("grid must be a list")

  if (scaled) scale <- list(shift=0, scale=1)  
  else        scale <- x$scale

  if (x$prior$priorK == "fixed"){  
    return(NMixPredCondDensJoint2.default(x=grid, icond=icond, scale=scale, K=x$K, w=as.numeric(t(x$w)), mu=as.numeric(t(x$mu)), Li=as.numeric(t(x$Li)), Krandom=FALSE))
  }else{
    return(NMixPredCondDensJoint2.default(x=grid, icond=icond, scale=scale, K=x$K, w=x$w, mu=x$mu, Li=x$Li, Krandom=TRUE))    
  }  
}  
