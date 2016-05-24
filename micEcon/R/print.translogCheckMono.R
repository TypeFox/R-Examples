print.translogCheckMono <- function( x, ... ) {

   if( length( x$increasing ) != ncol( x$exog ) ) {
      stop( "internal error: the length of 'x$increasing' is not equal to",
         " the number of columns of 'x$exog'" )
   }

   cat( "This translog function is" )
   if( sum( x$increasing ) > 0 ) {
      if( x$strict ) {
         cat( " strictly" )
      }
      cat( " monotonically increasing in " )
      cat( colnames( x$exog )[ x$increasing ], sep = ", " )
      if( sum( !x$increasing ) > 0 ) {
         cat( " and" )
      }
   }
   if( sum( !x$increasing ) > 0 ){
      if( x$strict ) {
         cat( " strictly" )
      }
      cat( " monotonically decreasing in " )
      cat( colnames( x$exog )[ !x$increasing ], sep = ", " )
   }
   cat( " at", sum( x$obs, na.rm = TRUE ) )
   cat( " out of", sum( !is.na( x$obs ) ), "observations (" )
   cat( round( 100 * sum( x$obs, na.rm = TRUE ) / sum( !is.na( x$obs ) ),
      digits = 1 ) )
   cat( "%)\n" )
   invisible( x )
}
