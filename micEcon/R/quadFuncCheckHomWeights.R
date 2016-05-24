.quadFuncCheckHomWeights <- function( homWeights, xNames ) {

   if( !is.null( homWeights ) ) {
      if( is.null( names( homWeights ) ) ) {
         stop( "the elements of argument 'homWeights' must have names" )
      }
      if( !all( names( homWeights ) %in% xNames ) ) {
         stop( "all names in argument 'homWeights' must be in argument 'xNames'" )
      }
      if( abs( sum( homWeights ) - 1 ) > .Machine$double.eps ^ 0.5 ) {
         stop( "the sum of the elements in argument 'homWeights' must be 1" )
      }
   }
   return( NULL )
}
