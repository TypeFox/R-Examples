print.summary.mvProbit <- function( x, digits = 4, ... ) {

   cat( "\n" )
   cat( "Call:\n" )
   cat( paste( deparse( x$call ), sep = "\n", collapse = "\n" ) )
   cat( "\n\n" )
   cat( "Coefficients:\n" )
   printCoefmat( coef( x ), digits = digits )
   cat( "\n" )
   cat( maximType( x ), ", ", nIter( x ), " iterations\n", sep = "" )
   cat( "Return code ", returnCode( x ), ": ", returnMessage( x ),
      "\n", sep = "" )
   cat( "Log-likelihood:", x$loglik, "on", sum( activePar( x ) ), "Df\n" )
   cat( "\n" )
   invisible( x )
}
