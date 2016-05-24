##
##  PURPOSE:   Computation of the predictive conditional univariate densities given one margin
##              * method for objects of class GLMM_MCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   24/07/2009
##
##  FUNCTIONS: NMixPredCondDensMarg.GLMM_MCMC (24/07/2009)
##             
## ======================================================================

## *************************************************************
## NMixPredCondDensMarg.GLMM_MCMC
## *************************************************************
NMixPredCondDensMarg.GLMM_MCMC <- function(x, icond, prob, grid, lgrid=50, scaled=FALSE, ...)
{
  if (missing(icond)) stop("icond must be given")
  if (x$dimb == 1) stop("not applicable for univariate mixtures")  
  
  if (missing(grid)){
    grid <- list()
    if (scaled){
      for (i in 1:x$dimb){
        rangeGrid <- 0 + c(-3.5, 3.5)*1
        grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
      }
    }else{
      for (i in 1:x$dimb){
        rangeGrid <- x$summ.b.Mean["Median", i] + c(-3.5, 3.5)*x$summ.b.SDCorr["Median", (i-1)*(2*x$dimb - i + 2)/2 + 1]
        grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
      }
    }
    names(grid) <- paste("b", 1:x$dimb, sep="")    
  }

  if (!is.list(grid)) stop("grid must be a list")

  if (scaled) scale <- list(shift=0, scale=1)  
  else        scale <- x$scale.b

  if (x$prior.b$priorK == "fixed"){  
    return(NMixPredCondDensMarg.default(x=grid, icond=icond, prob=prob, scale=scale, K=x$K_b, w=as.numeric(t(x$w_b)), mu=as.numeric(t(x$mu_b)), Li=as.numeric(t(x$Li_b)), Krandom=FALSE))
  }else{
    return(NMixPredCondDensMarg.default(x=grid, icond=icond, prob=prob, scale=scale, K=x$K_b, w=x$w_b, mu=x$mu_b, Li=x$Li_b, Krandom=TRUE))    
  }  
}  
