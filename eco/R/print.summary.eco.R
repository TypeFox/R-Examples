print.summary.eco <- function(x, digits=max(3, getOption("digits")-3), ...) {
    cat("\nCall: ") 
    cat(paste(deparse(x$call), sep="\n", collapse="\n"))

        cat("\n")
    if (!is.null(x$param.table)) {
           cat("\nParameter Estimates:\n")
           print(x$param.table, digits=digits, na.print="NA",...)
        }
 

   cat("\n*** Insample Predictions ***\n")
   cat("\nUnweighted:\n")
   print(x$agg.table, digits=digits, na.print="NA",...)
  
   if (!is.null(x$agg.wtable)) {
   cat("\nWeighted:\n")
   print(x$agg.wtable, digits=digits, na.print="NA",...)
  }
   
        cat("\nNumber of Units:", x$n.obs)
        cat("\nNumber of Monte Carlo Draws:", x$n.draws)
   
  
        if (!is.null(x$W1.table)) {
           cat("\n\nUnit-level Estimates of W1:\n")
           print(x$W1.table, digits=digits, na.print="NA",...)
           cat("\n\nUnit-level Estimates of W2:\n")
           print(x$W2.table, digits=digits, na.print="NA",...)
        }

      cat("\n")
      invisible(x)
}
