print.summary.aidsEst <- function( x, ... ) {
   cat( "\nDemand analysis with the Almost Ideal " )
   cat( "Demand System (AIDS)\n" )
   cat( "Estimation Method: " )
   cat( .aidsEstMethod( x$method, x$priceIndex ) )
   cat( "Estimated Coefficients:\n" )
   printCoefmat( x$coef$stat )

   cat( "R-squared Values of expenditure shares:\n" )
   print( x$r2 )

   cat( "R-squared Values of quantities:\n" )
   print( x$r2q )

   invisible( x )
}
