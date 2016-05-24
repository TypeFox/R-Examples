print.selection <- function( x,
      digits = max(3, getOption("digits") - 3), ... ) {
   cat( "\nCall:\n", deparse( x$call ), "\n\n" )
   cat( "Coefficients:\n" )
   print( coef( x ), digits = digits )
   cat( "\n" )
   invisible( x )
}