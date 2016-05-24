##
##  PURPOSE:   MCMC for a normal mixture model
##             * basic posterior summary for mixture components (after re-labeling to achieve
##               at least some identifiability) in a model with fixed number of components
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTION:  NMixSummComp.NMixMCMC
##              01/04/2015:  code allowing for a factor covariate on mixture weights
##                           added
##
## ======================================================================

## *************************************************************
## NMixSummComp.NMixMCMC
## *************************************************************
NMixSummComp.NMixMCMC <- function(x)
{
  #if (class(x) != "NMixMCMC") stop("x must be of class NMixMCMC")
  if (x$prior$priorK != "fixed") stop("number of mixture components was not fixed")

  if (x$dim == 1){
    for (k in 1:x$prior$Kmax){
      muk  <- x$poster.mean.mu[k,1] * x$scale$scale + x$scale$shift
      Vark <- x$scale$scale^2 * as.numeric(x$poster.mean.Sigma[[k]])
      SDk  <- sqrt(Vark)

      cat("\nComponent ", k, "\n", sep="")
      if (x$nx_w == 1){
        cat("    Weight:         ", x$poster.mean.w[k], "\n")
      }else{
        w2print <- x$poster.mean.w[seq(0, (x$nx_w - 1) * x$prior$Kmax, by = x$prior$Kmax) + k]
        names(w2print) <- x$lx_w
        cat("    Weights:\n")
        print(w2print)
        cat("\n")
      }    
      cat("    Mean:           ", muk, "\n")
      cat("    Variance:       ", Vark, "\n")
      cat("    Std. deviation: ", SDk, "\n")      
      cat("\n---------------------------------------------\n")
    }  
  }else{
    for (k in 1:x$prior$Kmax){
      muk  <- x$poster.mean.mu[k,] * x$scale$scale + x$scale$shift
      Vark <- diag(x$scale$scale) %*% x$poster.mean.Sigma[[k]] %*% diag(x$scale$scale)
      SDk  <- diag(sqrt(diag(Vark)))
      iSDk <- diag(1/diag(SDk))
      Cork <- iSDk %*% Vark %*% iSDk

      rownames(Vark) <- colnames(Vark) <- rownames(Cork) <- colnames(Cork) <- names(muk)
        
      cat("\nComponent ", k, "\n", sep="")
      if (x$nx_w == 1){      
        cat("    Weight: ", x$poster.mean.w[k], "\n")
      }else{
        w2print <- x$poster.mean.w[seq(0, (x$nx_w - 1) * x$prior$Kmax, by = x$prior$Kmax) + k]
        names(w2print) <- x$lx_w          
        cat("    Weights:\n")
        print(w2print)
        cat("\n")
      }    
      cat("    Mean:   ", muk, "\n")
      cat("\n")
      cat("    Covariance matrix:\n")
      print(Vark)
      cat("    Standard deviations: ", diag(SDk))
      cat("\n\n")
      cat("    Correlation matrix:\n")
      print(Cork)    
      cat("\n---------------------------------------------\n")
    }
  }  

  return(invisible(x))  
}
