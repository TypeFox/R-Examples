.snqProfitQuantNames <- function( quant, nNetput ) {
   if( !is.null( names( quant ) ) ) {
      quantNames <- names( quant )
   } else {
      quantNames <- paste( "q", 1:nNetput, sep = "" )
   }
   return( quantNames )
}

.snqProfitPriceNames <- function( prices, nNetput ) {
   if( !is.null( names( prices ) ) ) {
      priceNames <- names( prices )
   } else {
      priceNames <- paste( "p", 1:nNetput, sep = "" )
   }
   return( priceNames )
}
