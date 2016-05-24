print.summary.censReg <- function( x, logSigma = TRUE, digits = 4, ... ) {

   cat( "\n" )
   cat( "Call:\n" )
   cat( paste( deparse( x$call ), sep = "\n", collapse = "\n" ) )
   cat( "\n\n" )
   cat( "Observations:\n" )
   print( x$nObs )
   cat( "\n" )
   cat( "Coefficients:\n" )
   printCoefmat( coef( x, logSigma = logSigma ), digits = digits )
   cat( "\n" )
   cat( maximType( x ), ", ", nIter( x ), " iterations\n", sep = "" )
   cat( "Return code ", returnCode( x ), ": ", returnMessage( x ),
      "\n", sep = "" )
   cat( "Log-likelihood:", x$loglik, "on", sum( activePar( x ) ), "Df\n" )
   cat( "\n" )
   invisible( x )
}
