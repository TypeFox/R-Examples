print.aidsEst <- function( x, ... ) {
   cat( "\nDemand analysis with the Almost Ideal " )
   cat( "Demand System (AIDS)\n" )
   cat( "Estimation Method: " )
   cat( .aidsEstMethod( x$method, x$priceIndex ) )
   cat( "Coefficients:\n" )
   print( coef( x ) )
   invisible( x )
}
