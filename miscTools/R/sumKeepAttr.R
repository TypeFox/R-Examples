## return the sum of an array while keeping its attributes
sumKeepAttr <- function( x, keepNames = FALSE, na.rm = FALSE ) {
   xAttr <- attributes( x )
   if( !keepNames ) {
      xAttr$names <- NULL
   }
   x <- sum( x, na.rm = na.rm )
   mostattributes( x ) <- xAttr
   return( x )
}
