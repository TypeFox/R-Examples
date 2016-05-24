##
##  PURPOSE:   MCMC for a GLMM with normal mixture in the random effects distribution
##             * basic posterior summary for mixture components (after re-labeling to achieve
##               at least some identifiability) in a model with fixed number of components
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   17/06/2009
##
##  FUNCTION:  NMixSummComp.GLMM_MCMC
##             
##
## ======================================================================

## *************************************************************
## NMixSummComp.GLMM_MCMC
## *************************************************************
NMixSummComp.GLMM_MCMC <- function(x)
{
  #if (class(x) != "GLMM_MCMC") stop("x must be of class GLMM_MCMC")
  if (x$prior.b$priorK != "fixed") stop("number of mixture components was not fixed")

  if (x$dimb == 1){
    for (k in 1:x$prior.b$Kmax){
      muk  <- x$poster.mean.mu_b[k,1] * x$scale.b$scale + x$scale.b$shift
      Vark <- x$scale.b$scale^2 * as.numeric(x$poster.mean.Sigma_b[[k]])
      SDk  <- sqrt(Vark)
        
      cat("\nComponent ", k, "\n", sep="")
      cat("    Weight:         ", x$poster.mean.w_b[k], "\n")
      cat("    Mean:           ", muk, "\n")
      cat("    Variance:       ", Vark, "\n")
      cat("    Std. deviation: ", SDk, "\n")      
      cat("\n---------------------------------------------\n")
    }    
  }else{  
    for (k in 1:x$prior.b$Kmax){
      muk  <- x$poster.mean.mu_b[k,] * x$scale.b$scale + x$scale.b$shift
      Vark <- diag(x$scale.b$scale) %*% x$poster.mean.Sigma_b[[k]] %*% diag(x$scale.b$scale)
      SDk  <- diag(sqrt(diag(Vark)))
      iSDk <- diag(1/diag(SDk))
      Cork <- iSDk %*% Vark %*% iSDk

      rownames(Vark) <- colnames(Vark) <- rownames(Cork) <- colnames(Cork) <- names(muk)
        
      cat("\nComponent ", k, "\n", sep="")
      cat("    Weight: ", x$poster.mean.w_b[k], "\n")
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
