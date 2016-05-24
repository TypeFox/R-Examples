summary.mvProbit <- function( object, ... ) {

   result <- NextMethod( "summary", object )

   result$call <- object$call
   result$start <- object$start
   result$nDep <- object$nDep
   result$nReg <- object$nReg
   result$nObs <- object$nObs

   class( result ) <- c( "summary.mvProbit", class( result ) )

   return( result )
}
