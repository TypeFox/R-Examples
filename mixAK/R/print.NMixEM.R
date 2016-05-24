##
##  PURPOSE:   EM algorithm to compute ML estimates in a normal mixture
##             * print method for the resulting object
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:                  31/12/2009
##
##  FUNCTIONS:  print.NMixEM
##
## ======================================================================

## *************************************************************
## print.NMixEM
## *************************************************************
print.NMixEM <- function(x, ...)
{
  cat("\n")
  cat(paste("     ", x$K, " component normal mixture estimated using EM algorithm\n", sep=""))
  cat("     =======================================================\n\n")
  
  if (x$dim == 1){
    SD <- sqrt(x$Sigma)
    cat("Component variance:       ", as.numeric(x$Sigma), "\n")
    cat("Component std. deviation: ", SD, "\n")      
    cat("\n---------------------------------------------\n")
    
    for (k in 1:x$K){        
      cat("\nComponent ", k, "\n", sep="")
      cat("    Weight:         ", x$weight[k], "\n")
      cat("    Mean:           ", x$mean[k], "\n")
      cat("\n---------------------------------------------\n")
    }    
  }else{
    SD <- diag(sqrt(diag(x$Sigma)))
    iSD <- diag(1/diag(SD))
    Cor <- iSD %*% x$Sigma %*% iSD
    rownames(x$Sigma) <- colnames(x$Sigma) <- rownames(Cor) <- colnames(Cor) <- colnames(x$mean)
    
    cat("Component covariance matrix:\n")
    print(x$Sigma)
    cat("Component standard deviations: ", diag(SD))
    cat("\n\n")
    cat("Component correlation matrix:\n")
    print(Cor)    
    cat("\n---------------------------------------------\n")
    
    for (k in 1:x$K){        
      cat("\nComponent ", k, "\n", sep="")
      cat("    Weight: ", x$weight[k], "\n")
      cat("    Mean:   ", x$mean[k,], "\n")
      cat("\n---------------------------------------------\n")
    }
  }  

  cat("\n")
  cat("Log-likelihood: ", x$loglik, ",  AIC: ", x$aic, ", BIC: ", x$bic, "\n", sep="")
  cat("EM iterations:  ", x$iter, "\n")
  cat("\n")
  
  return(invisible(x))  
}
