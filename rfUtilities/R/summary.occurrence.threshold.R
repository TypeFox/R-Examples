#' @title Summarizing occurrence.threshold
#' @description Summarize occurrence.threshold 
#'
#' @param object  Object of occurrence.threshold
#' @param ... Ignored
#'
#' @method summary occurrence.threshold
#'
#' @export
summary.occurrence.threshold <- function(object, ...) { 
  cat("Evaluation statistic:", object$statistic, "\n")
    cat("\n")  
    cat("Moments for thresholds:", "\n")
    cat("\n")  
    print( summary(object$thresholds) )
    cat("\n")  	
  if(object$statistic == "delta.ss") { 	
    cat("Probability threshold with minimum delta sensitivity/specificity = ", 
      names(object$thresholds)[min(which(object$thresholds == min(object$thresholds)))], "\n")
	  } 
	else if (object$statistic == "sum.ss") {
    cat("Probability threshold with maximum cummlative sensitivity/specificity = ", 
      names(object$thresholds)[min(which(object$thresholds == max(object$thresholds)))], "\n")
      }
    else if (object$statistic == "kappa") {
    cat("Probability threshold with maximum kappa = ", 
      names(object$thresholds)[min(which(object$thresholds == max(object$thresholds)))], "\n")	  
  }	  
} 
