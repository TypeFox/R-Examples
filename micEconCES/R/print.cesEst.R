print.cesEst <- function( x, digits = max( 3, getOption( "digits" ) - 3 ),
      ... ) {

   cat( "Estimated CES function\n\n" )

   cat( "Call:\n" )
   print( x$call )
   cat( "\n" )

   cat( "Coefficients:\n" )
   print( coef( x ), digits = digits )
   cat( "\n" )

   if( !is.null( x$ela ) ) {
      if( length( x$ela ) == 1 ) {
         cat( "Elasticity of Substitution:", 
            format( x$ela, digits = digits ), "\n" )
      } else {
         cat( "Elasticities of Substitution:\n" )
         print( x$ela, digits = digits )
         cat( "AU = Allen-Uzawa (partial) elasticity of substitution\n" )
         cat( "HM = Hicks-McFadden (direct) elasticity of substitution\n" )
      }
      cat( "\n" )
   }

   invisible( x )
}
