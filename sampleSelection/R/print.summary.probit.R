print.summary.probit <- function( x, ... ) {
   cat("--------------------------------------------\n")
   cat("Probit binary choice model/Maximum Likelihood estimation\n")
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", logLik(x), "\n")
      cat("Model: Y == '", x$levels[2], "' in contrary to '", x$levels[1], "'\n", sep="")
      cat(x$nObs, " observations (", x$N0, " 'negative' and ", x$N1, " 'positive') and ",
          x$NActivePar, " free parameters (df = ",
          x$df.residual, ")\n", sep="")
      cat("Estimates:\n")
      printCoefmat( x$estimate, ... )
   }
   cat("Significance test:\n")
   cat("chi2(", x$LRT$df, ") = ", x$LRT$LRT, " (p=", x$LRT$pchi2, ")\n", sep="")
   cat("--------------------------------------------\n")
}
