print.coef.selection <- function( x, prefix = TRUE,
      digits = max(3, getOption("digits") - 3), ... ) {

   addToCoefNames <- function( coefValues, prefix, index ) {
      if( !is.null( index ) ) {
         names( coefValues )[ index ] <-
            paste( prefix, names( coefValues )[ index ], sep = "" )
      }
      return( coefValues )
   }
   if( attributes( x )$part == "full" && prefix ){
      if( attributes( x )$tobitType == 2) {
         x <- addToCoefNames( x, "S:", attributes( x )$index$betaS )
         x <- addToCoefNames( x, "O:", attributes( x )$index$betaO )
      } else if( attributes( x )$tobitType == 5) {
         x <- addToCoefNames( x, "S:",  attributes( x )$index$betaS )
         x <- addToCoefNames( x, "O1:", attributes( x )$index$betaO1 )
         x <- addToCoefNames( x, "O2:", attributes( x )$index$betaO2 )
      }
   }

   print.default( format( x, digits = digits ),
      print.gap = 2, quote = FALSE )
   invisible( x )
}
