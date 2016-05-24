print.translogEst <- function( x, ... ) {

   cat( "Estimated Translog function with", x$nObs, "observations.\n" )
   print( coef( x ) )

   invisible( x )
}
