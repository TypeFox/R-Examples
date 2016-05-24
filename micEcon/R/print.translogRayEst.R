print.translogRayEst <- function( x, ... ) {

   cat( "Estimated Translog ray function with ",
      length( x$yNames ), " dependent variables,\n",
      length( x$xNames ), " independent variables, and ",
      x$nObs, " observations.\n", sep = "" )
   print( coef( x ) )

   invisible( x )
}
