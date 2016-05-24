print.frontier <- function( x, digits = NULL, ... ) {

   if( is.null( digits ) ) {
      digits <- max( 3, getOption( "digits" ) - 3 )
   }
   cat( "\nCall:\n" )
   cat( deparse( x$call ) )
   cat( "\n\n" )
   cat( "Maximum likelihood estimates\n" )
   print.default( format( coef( x ), digits = digits ), print.gap = 2,
      quote = FALSE )
   invisible( x )
}