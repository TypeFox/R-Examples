print.maxLik <- function( x, ... ) {
   cat("Maximum Likelihood estimation\n")
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", x$maximum )
      cat( " (", sum( activePar( x ) ), " free parameter(s))\n", sep = "" )
      cat("Estimate(s):", x$estimate, "\n" )
   }
}
