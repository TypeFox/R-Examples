#' @title Print significance
#' @description print method for class "significance"
#' @param x    Object of class significance
#' @param ...  Ignored
#'
#' @method print significance
#'
#' @export
print.significance <- function(x, ...) {
cat("Number of permutations: ", x$nPerm, "\n")
  cat("p-value: ", round(x$pValue,3), "\n")
  if( x$pValue <= x$pValueThreshold ) accept=TRUE else accept=FALSE 
    if (accept == TRUE) accept <- paste("Model signifiant at p = ", round(x$pValue,3), sep="" ) 
	if (accept == FALSE) accept <- paste("Model not signifiant at p = ", round(x$pValue,3), sep="" )
    cat(accept, "\n")
	  if(x$rf.type == "classification") {
	    cat("\t", "Model OOB error: ", x$test.OOB, "\n")
	    cat("\t", "Random OOB error: ", stats::median(x$RandOOB), "\n")
	    cat("\t", "min random global error:", min(x$RandOOB), "\n")
	    cat("\t", "max random global error: ",  max(x$RandOOB), "\n")
	    cat("\t", "min random within class error:", min(x$RandMaxError), "\n")
	    cat("\t", "max random within class error: ", min(x$RandMaxError), "\n")
	  } else if(x$rf.type == "regression") {
	    cat("\t", "Model R-square: ", x$Rsquare, "\n")
	    cat("\t", "Random R-square: ", stats::median(x$RandRsquare), "\n")
	    cat("\t", "Random R-square variance: ", stats::var(x$RandRsquare), "\n")		
	  }
}	
