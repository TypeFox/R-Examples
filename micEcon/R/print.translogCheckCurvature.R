print.translogCheckCurvature <- function( x, ... ) {

   cat( "This translog function is " )
   if( x$quasi ) {
      cat( "quasi" )
   }
   if( x$convex ) {
      cat( "convex" )
   } else {
      cat( "concave" )
   }
   cat( " at", sum( x$obs, na.rm = TRUE ) )
   cat( " out of", sum( !is.na( x$obs ) ), "observations (" )
   cat( round( 100 * sum( x$obs, na.rm = TRUE ) / sum( !is.na( x$obs ) ),
      digits = 1 ) )
   cat( "%)\n" )
   invisible( x )
}
