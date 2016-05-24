##
##  PURPOSE:   Computation of the conditional univariate densities given one margin
##             (plug-in version)
##              * method for objects of class NMixMCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTIONS: NMixPlugCondDensMarg.NMixMCMC (28/05/2009)
##             
## ======================================================================

## *************************************************************
## NMixPlugCondDensMarg.NMixMCMC
## *************************************************************
NMixPlugCondDensMarg.NMixMCMC <- function(x, icond, grid, lgrid=50, scaled=FALSE, ...)
{
  if (missing(icond)) stop("icond must be given")
  if (x$dim == 1) stop("not applicable for univariate mixtures")  
  if (x$prior$priorK != "fixed") stop("only implemented for models with fixed number of components")

  if (x$nx_w > 1) stop("This function has not (yet) been implemented if a factor covariate on mixture weights is present.")
  
  if (missing(grid)){
    grid <- list()
    if (scaled){
      if (x$dim == 1){
        rangeGrid <- x$summ.z.Mean["Median"] + c(-3.5, 3.5)*x$summ.z.SDCorr["Median"]
        grid[[1]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)      
      }else{     
        for (i in 1:x$dim){
          rangeGrid <- x$summ.z.Mean["Median", i] + c(-3.5, 3.5)*x$summ.z.SDCorr["Median", (i-1)*(2*x$dim - i + 2)/2 + 1]
          grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
        }
      }  
    }else{
      if (x$dim == 1){
        rangeGrid <- x$summ.y.Mean["Median"] + c(-3.5, 3.5)*x$summ.y.SDCorr["Median"]
        grid[[1]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)      
      }else{     
        for (i in 1:x$dim){
          rangeGrid <- x$summ.y.Mean["Median", i] + c(-3.5, 3.5)*x$summ.y.SDCorr["Median", (i-1)*(2*x$dim - i + 2)/2 + 1]
          grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
        }
      }  
    }
    names(grid) <- paste("x", 1:x$dim, sep="")    
  }

  if (x$dim == 1) if (is.numeric(grid)) grid <- list(x1=grid)
  if (!is.list(grid)) stop("grid must be a list")

  if (scaled) scale <- list(shift=0, scale=1)  
  else        scale <- x$scale
  
  return(NMixPlugCondDensMarg.default(x=grid, icond=icond, scale=scale, w=x$poster.mean.w, mu=x$poster.mean.mu, Sigma=x$poster.mean.Sigma))
}  
