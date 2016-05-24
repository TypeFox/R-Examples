#' @title Print occurrence.threshold
#' @description Print method for occurrence.threshold objects
#'    
#' @param x    Object of class occurrence.threshold
#' @param ...  Ignored
#'
#' @method print occurrence.threshold
#'
#' @export
print.occurrence.threshold <- function(x, ...) {
  cat("Evaluation statistic:", x$statistic, "\n")
    cat("\n")  
    cat("Moments for thresholds:", "\n")
    cat("\n")  
    print( summary(x$thresholds) )
    cat("\n")  	
  if(x$statistic == "delta.ss") { 	
    cat("Probability threshold with minimum delta sensitivity/specificity = ", 
      names(x$thresholds)[min(which(x$thresholds == min(x$thresholds)))], "\n")
	  } 
	else if (x$statistic == "sum.ss") {
    cat("Probability threshold with maximum cummlative sensitivity/specificity = ", 
      names(x$thresholds)[min(which(x$thresholds == max(x$thresholds)))], "\n")
      }
    else if (x$statistic == "kappa") {
    cat("Probability threshold with maximum kappa = ", 
      names(x$thresholds)[min(which(x$thresholds == max(x$thresholds)))], "\n")	  
  }	  
} 
