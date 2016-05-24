colMedians <- function( x, na.rm = FALSE ) {

   if( is.data.frame( x ) ) {
      x <- as.matrix( x )
   }
   if( !is.array( x ) ) {
      stop( "argument 'x' must be a data.frame, matrix, or array" )
   }
   if( !is.numeric( x ) ) {
      stop( "argument 'x' must be numeric" )
   }

   result <- array( NA, dim = dim( x )[-1] )
   dimnames( result ) <- dimnames( x )[-1]

   for( i in 1:dim( x )[ 2 ] ) {
      if( length( dim( x ) ) == 2 ) {
         result[ i ] <- median( x[ , i ], na.rm = na.rm )
      } else {
         result[ slice.index( result, 1 ) == i ] <-
            colMedians( array( x[ slice.index( x, 2 ) == i ],
               dim = dim( x )[ -2 ] ), na.rm = na.rm )
      }
   }

   return( result )
}

rowMedians <- function( x, na.rm = FALSE ) {
   if( is.null( dim( x ) ) || length( dim( x ) ) != 2 ) {
      stop( "argument 'x' must be a matrix or a data.frame" )
   }
   return( colMedians( t( x ), na.rm = na.rm ) )
}
