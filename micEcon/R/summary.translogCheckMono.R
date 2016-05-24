summary.translogCheckMono <- function( object, ... ) {

   class( object ) <- "summary.translogCheckMono"
   return( object )
}

print.summary.translogCheckMono <- function( x, ... ) {

   print.translogCheckMono( x )

   cat( "The monotonicity condition for the exogenous variable\n" )
   for( i in 1:ncol( x$exog ) ) {
      cat( "- '", colnames( x$exog )[ i ], "'", sep = "" )
      cat( " is fulfilled at", sum( x$exog[ , i ], na.rm = TRUE ) )
      cat( " out of", sum( !is.na( x$exog[ , i ] ) ), "observations (" )
      cat( round( 100 * sum( x$exog[ , i ], na.rm = TRUE ) /
         sum( !is.na( x$exog[ , i ] ) ), digits = 1 ) )
      cat( "%)\n" )
   }
   invisible( x )
}
