##
##  PURPOSE:   Computation of the predictive conditional univariate densities given one margin
##              * method for objects of class NMixMCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   31/05/2009
##
##  FUNCTIONS: NMixPredCondDensMarg.NMixMCMC (31/05/2009)
##             
## ======================================================================

## *************************************************************
## NMixPredCondDensMarg.NMixMCMC
## *************************************************************
NMixPredCondDensMarg.NMixMCMC <- function(x, icond, prob, grid, lgrid=50, scaled=FALSE, ...)
{
  if (missing(icond)) stop("icond must be given")
  if (x$dim == 1) stop("not applicable for univariate mixtures")  

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
    return(NMixPredCondDensMarg.default(x=grid, icond=icond, prob=prob, scale=scale, K=x$K, w=as.numeric(t(x$w)), mu=as.numeric(t(x$mu)), Li=as.numeric(t(x$Li)), Krandom=FALSE))
  }else{
    return(NMixPredCondDensMarg.default(x=grid, icond=icond, prob=prob, scale=scale, K=x$K, w=x$w, mu=x$mu, Li=x$Li, Krandom=TRUE))    
  }  
}  
