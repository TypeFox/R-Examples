##
##  PURPOSE:   Generalized linear mixed model with possibly several response variables
##             and normal mixtures in the distribution of continuous response and random effects
##             * print method for the resulting object
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    18/07/2009
##
##  FUNCTIONS:  print.GLMM_MCMC
##
## ================================================================================================

## *************************************************************
## print.GLMM_MCMC
## *************************************************************
print.GLMM_MCMC <- function(x, ...)
{
  R <- sum(x$R)
  
  cat("\n")
  cat(paste("     Generalized linear mixed model for ", R, " responses estimated using MCMC\n", sep=""))
  cat("     ====================================================================\n")

  cat("\nDeviance posterior summary statistics:")
  cat("\n-----------------------------------------------\n")
  pr.Dev <- x$summ.Deviance[, "Deviance"]
  names(pr.Dev) <- rownames(x$summ.Deviance)
  print(pr.Dev, ...)
      
  if (x$lalpha){
    cat("\nPosterior summary statistics for fixed effects:")
    cat("\n-----------------------------------------------\n")
    print(x$summ.alpha, ...)
  }

  if (x$dimb){
    if (x$prior.b$priorK == "fixed"){
      cat(paste("\nDistribution of random effects is a normal mixture with ", x$prior.b$Kmax, " components", sep=""))
      cat("\n---------------------------------------------------------------------")
    }else{
      cat(paste("\nDistribution of random effects is a normal mixture with at most ", x$prior.b$Kmax, " components", sep=""))
      cat("\n------------------------------------------------------------------------------------------------")
      cat("\nPosterior distribution of K (number of mixture components):")
      cat("\n-----------------------------------------------------------")
      print(x$propK_b, ...)    
    }
    cat("\nPosterior summary statistics for moments of mixture for random effects:")
    cat("\n-----------------------------------------------------------------------")
    if (x$dimb == 1){
      cat("\nMean:\n")
      print(x$summ.b.Mean, ...)

      cat("\nStandard deviation:\n")
      print(x$summ.b.SDCorr, ...)
    }else{
      cat("\nMeans:\n")
      print(x$summ.b.Mean, ...)

      cat("\nStandard deviations and correlations:\n")
      print(x$summ.b.SDCorr, ...)
    }      
  }

  if (x$R["Rc"]){
    cat("\nPosterior summary statistics for standard deviations")
    cat("\nof residuals of continuous responses:")
    cat("\n----------------------------------------------------\n")
    print(x$summ.sigma_eps, ...)
  }  

  return(invisible(x))  
}
