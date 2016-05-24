##
##  PURPOSE:   Computation of the predictive marginal (univariate) cumulative distribution functions
##              * method for objects of class GLMM_MCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   24/07/2009
##
##  FUNCTIONS: NMixPredCDFMarg.GLMM_MCMC (24/07/2009)
##             
## ======================================================================

## *************************************************************
## NMixPredCDFMarg.GLMM_MCMC
## *************************************************************
NMixPredCDFMarg.GLMM_MCMC <- function(x, grid, lgrid=500, scaled=FALSE, ...)
{
  if (missing(grid)){
    grid <- list()
    if (scaled){
      if (x$dimb == 1){
        rangeGrid <- 0 + c(-3.5, 3.5)*1
        grid[[1]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)      
      }else{     
        for (i in 1:x$dimb){
          rangeGrid <- 0 + c(-3.5, 3.5)*1
          grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
        }
      }  
    }else{
      if (x$dimb == 1){
        rangeGrid <- x$summ.b.Mean["Median"] + c(-3.5, 3.5)*x$summ.b.SDCorr["Median"]
        grid[[1]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)      
      }else{     
        for (i in 1:x$dimb){
          rangeGrid <- x$summ.b.Mean["Median", i] + c(-3.5, 3.5)*x$summ.b.SDCorr["Median", (i-1)*(2*x$dimb - i + 2)/2 + 1]
          grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
        }
      }  
    }
    names(grid) <- paste("b", 1:x$dimb, sep="")    
  }

  if (x$dimb == 1) if (is.numeric(grid)) grid <- list(x1=grid)
  if (!is.list(grid)) stop("grid must be a list")

  if (scaled) scale <- list(shift=0, scale=1)  
  else        scale <- x$scale.b
  
  if (x$prior.b$priorK == "fixed"){
    return(NMixPredCDFMarg.default(x=grid, scale=scale, K=x$K_b, w=as.numeric(t(x$w_b)), mu=as.numeric(t(x$mu_b)), Li=as.numeric(t(x$Li_b)), Krandom=FALSE))
  }else{
    return(NMixPredCDFMarg.default(x=grid, scale=scale, K=x$K_b, w=x$w_b, mu=x$mu_b, Li=x$Li_b, Krandom=TRUE))
  }    
}  
