summary.censReg <- function( object, ... ) {

   result <- NextMethod( summary, object )

   result$call <- object$call
   result$nObs <- object$nObs

   class( result ) <- c( "summary.censReg", class( result ) )

   return( result )
}
