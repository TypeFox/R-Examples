print.mvProbit <- function( x, digits = 4, ... ) {

   cat( "\n" )
   cat( "Call:\n" )
   cat( paste( deparse( x$call ), sep = "\n", collapse = "\n" ) )
   cat( "\n\n" )
   cat( "Coefficients:\n" )
   print( coef( x ), digits = digits )
   cat( "\n" )
   invisible( x )
}
