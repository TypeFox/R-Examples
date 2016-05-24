##
##  PURPOSE:   Computation of the predictive marginal (univariate) cumulative distribution functions
##              * method for objects of class NMixMCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/06/2009
##             02/04/2015:  code allowing for a factor covariate
##                          on mixture weights implemented
##
##  FUNCTIONS: NMixPredCDFMarg.NMixMCMC (09/06/2009)
##             
## ======================================================================

## *************************************************************
## NMixPredCDFMarg.NMixMCMC
## *************************************************************
NMixPredCDFMarg.NMixMCMC <- function(x, grid, lgrid=500, scaled=FALSE, ...)
{
  if (scaled) scale <- list(shift=0, scale=1)  
  else        scale <- x$scale

  if (x$nx_w == 1){   
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
    }

    if (x$dim == 1) if (is.numeric(grid)) grid <- list(x1=grid)
    if (!is.list(grid)) stop("grid must be a list")
    names(grid) <- paste("x", 1:x$dim, sep="")    
    
    if (x$prior$priorK == "fixed"){
      cdfm <- NMixPredCDFMarg.default(x=grid, scale=scale, K=x$K, w=as.numeric(t(x$w)), mu=as.numeric(t(x$mu)), Li=as.numeric(t(x$Li)), Krandom=FALSE)
    }else{
      cdfm <- NMixPredCDFMarg.default(x=grid, scale=scale, K=x$K, w=x$w, mu=x$mu, Li=x$Li, Krandom=TRUE)
    }

    cdfm$nx_w <- 1
    cdfm$lx_w <- ""
    
    return(cdfm)
    
  }

  else{
    if (x$prior$priorK != "fixed") stop("only implemented for models with fixed number of components")
      
    if (missing(grid)){
      grid <- list()
      if (scaled){
        for (ixw in 1:x$nx_w){
          for (i in 1:x$dim){
            naamM <- paste("z.Mean.", i, "-", x$lx_w[ixw], sep = "")
            naamS <- paste("z.SD.", i, "-", x$lx_w[ixw], sep = "")
            rangeGrid <- x$summ.z.Mean["Median", naamM] + c(-3.5, 3.5) * x$summ.z.SDCorr["Median", naamS]
            grid[[(ixw - 1)*x$dim + i]] <- seq(rangeGrid[1], rangeGrid[2], length = lgrid)
          }    
        }    
      }else{
        for (ixw in 1:x$nx_w){
          for (i in 1:x$dim){
            naamM <- paste("y.Mean.", i, "-", x$lx_w[ixw], sep = "")
            naamS <- paste("y.SD.", i, "-", x$lx_w[ixw], sep = "")
            rangeGrid <- x$summ.y.Mean["Median", naamM] + c(-3.5, 3.5) * x$summ.y.SDCorr["Median", naamS]
            grid[[(ixw - 1)*x$dim + i]] <- seq(rangeGrid[1], rangeGrid[2], length = lgrid)
          }    
        }    
      }
    }    

    if (!is.list(grid)) stop("grid must be a list")
    names(grid) <- paste("x", rep(1:x$dim, x$nx_w), "-", rep(x$lx_w, each = x$dim), sep="")
    
    cdfm <- list(x = grid)
    for (ixw in 1:x$nx_w){     
      gridixw <- list()
      for (i in 1:x$dim){
        gridixw[[i]] <- grid[[(ixw - 1)*x$dim + i]]
      }
      names(gridixw) <- paste("x", 1:x$dim, sep = "")

      tmp <- NMixPredCDFMarg.default(x = gridixw, scale = scale, K = x$K, w = as.numeric(t(x$w[, (ixw - 1)*x$K[1] + 1:x$K[1]])), mu = as.numeric(t(x$mu)), Li = as.numeric(t(x$Li)), Krandom = FALSE)
      if (ixw == 1){
        cdfm$cdf <- tmp$cdf
      }else{
        for (i in 1:x$dim) cdfm$cdf[[(ixw - 1) * x$dim + i]] <- tmp$cdf[[i]]
      }    
    }
    names(cdfm$cdf) <- paste(rep(1:x$dim, x$nx_w), "-", rep(x$lx_w, each = x$dim), sep="")
    cdfm$nx_w <- x$nx_w
    cdfm$lx_w <- x$lx_w
    
    class(cdfm) <- "NMixPredCDFMarg"
    return(cdfm)           
  }    
}  


