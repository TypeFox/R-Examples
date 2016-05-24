print.synchrony <- function (x, digits=max(3L, getOption("digits") - 3L), ...) {
  # Community synchrony
  if (is.null(x$w.corrected) & !is.null(x$meancorr)) {
    cat("Community synchrony: ")
    cat(format(x$obs, digits=digits))
    
    cat("\nMean pairwise correlation: ")
    cat(format(x$meancorr, digits=digits))
    
    if (!is.null(x$pval)) {
      tail=paste0("one-tailed test [", x$alternative, "]")
      cat(paste("\nCommunity synchrony p-value (", tail, "): ", sep=""))
      cat(format(x$pval, digits=digits))
    }
  }
  
  # Mean correlation
  else if (is.null(x$w.corrected) & is.null(x$meancorr) & is.null(x$index)) {
    cat(paste("Mean", x$method, "correlation: ", sep=" "))
    cat(format(x$obs, digits=digits))
        
    if (!is.null(x$pval)) {
      tail=ifelse(x$alternative=="two.tailed", "two-tailed test", paste0("one-tailed test [", x$alternative, "]"))
      cat(paste("\nMean correlation p-value (", tail, "): ", sep=""))
      cat(format(x$pval, digits=digits))
    }
  }
  
  else if (!is.null(x$index)) {
    cat("Proportion of common peaks: ")
    cat(format(x$obs, digits=digits))
        
    if (!is.null(x$pval)) {
      cat("\nPeaks p-value (one-tailed test [greater]): ")
      cat(format(x$pval, digits=digits))
    }    
  }
  
  # Kendall's W
  else {
    cat("Kendall's W (uncorrected for ties): ")
    cat(format(x$w.uncorrected, digits=digits))
    
    cat("\nKendall's W (corrected for ties): ")
    cat(format(x$w.corrected, digits=digits))
    
    cat("\nSpearman's ranked correlation: ")
    cat(format(x$spearman.corr, digits=digits))
    
    if (!is.null(x$pval.rand)) {
      cat("\nKendall's W p-value (one-tailed test [greater]): ")
      cat(format(x$pval.rand, digits=digits))
    }
  }
}
