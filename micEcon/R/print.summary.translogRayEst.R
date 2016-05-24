print.summary.translogRayEst <- function( x, ... ) {

   cat( "Estimated Translog ray function with ",
      length( x$yNames ), " dependent variables,\n",
      length( x$xNames ), " independent variables, and ",
      x$nObs, " observations.\n", sep = "" )
   printCoefmat( x$coefTable )

   cat( "R-squared:", x$r2, "    Adjusted R-squared:", x$r2bar, "\n" )

   invisible( x )
}
