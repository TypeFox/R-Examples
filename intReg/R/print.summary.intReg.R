
print.summary.intReg <- function(x,
                                 digits=max(3, getOption("digits") - 3),
                                 ...) {
   cat("--------------------------------------------\n")
   cat("Interval regression\n")
   cat( "Maximum Likelihood estimation\n" )
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", logLik(x), "\n")
   }
   if(!is.null(x$estimate)) {
      cat( nObs(x), "observations, " )
      cat( sum(activePar(x)), "free parameters" )
      cat( " (df = ", x$param$df, ")\n", sep="")
      printCoefmat( x$estimate, signif.legend = TRUE, digits = digits )
   }
   cat("--------------------------------------------\n")
   invisible( x )
}
