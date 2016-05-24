print.summary.translogEst <- function( x, ... ) {

   cat( "Estimated Translog function with", x$nObs, "observations.\n" )
   printCoefmat( x$coefTable )

   cat( "R-squared:", x$r2, "    Adjusted R-squared:", x$r2bar, "\n" )

   invisible( x )
}
