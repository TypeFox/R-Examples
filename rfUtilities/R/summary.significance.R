#' @title Summarizing significance
#' @description Summarizing of a rf.significance object
#' @param object  Object of class significance
#' @param ... Ignored
#'
#' @method summary significance
#' @export
summary.significance <- function(object, ...) {
cat("Number of permutations: ", object$nPerm, "\n")
  cat("p-value: ", round(object$pValue,3), "\n")
  if( object$pValue <= object$pValueThreshold ) accept=TRUE else accept=FALSE 
    if (accept == TRUE) accept <- paste("Model signifiant at p = ", round(object$pValue,3), sep="" ) 
	if (accept == FALSE) accept <- paste("Model not signifiant at p = ", round(object$pValue,3), sep="" )
    cat(accept, "\n")
	  if(object$rf.type == "classification") {
	    cat("\t", "Model OOB error: ", object$test.OOB, "\n")
	    cat("\t", "Random OOB error: ", stats::median(object$RandOOB), "\n")
	    cat("\t", "min random global error:", min(object$RandOOB), "\n")
	    cat("\t", "max random global error: ",  max(object$RandOOB), "\n")
	    cat("\t", "min random within class error:", min(object$RandMaxError), "\n")
	    cat("\t", "max random within class error: ", min(object$RandMaxError), "\n")
	  } else if(object$rf.type == "regression") {
	    cat("\t", "Model R-square: ", object$Rsquare, "\n")
	    cat("\t", "Random R-square: ", stats::median(object$RandRsquare), "\n")
	    cat("\t", "Random R-square variance: ", stats::var(object$RandRsquare), "\n")		
	  }
}
	  