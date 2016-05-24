.aidsQuantNames <- function( shares, coef, nGoods ) {
   if( !is.null( names( shares ) ) ) {
      qNames <- paste( "q_", names( shares ), sep = "" )
   } else if( !is.null( rownames( coef$gamma ) ) ) {
      qNames <- paste( "q_", rownames( coef$gamma ), sep = "" )
   } else {
      qNames <- paste( "q", c( 1:nGoods ), sep = "" )
   }
   return( qNames )
}

.aidsPriceNames <- function( prices, coef, nGoods ) {
   if( !is.null( names( prices ) ) ) {
      priceNames <- names( prices )
   } else if( !is.null( colnames( coef$gamma ) ) ) {
      priceNames <- colnames( coef$gamma )
   } else {
      priceNames <- paste( "p", c( 1:nGoods ), sep = "" )
   }
   return( priceNames )
}
