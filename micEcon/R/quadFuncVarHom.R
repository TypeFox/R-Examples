.quadFuncVarHom <- function( data, xName, homWeights, deflator, 
      xSubtract = NULL ) {

   if( is.null( homWeights ) | ! xName %in% names( homWeights ) ) {
      result <- data[[ xName ]]
   } else {
      result <- data[[ xName ]] / deflator
      if( !is.null( xSubtract ) ) {
         result <- result - data[[ xSubtract ]] / deflator
      }
   }

   return( result )
}