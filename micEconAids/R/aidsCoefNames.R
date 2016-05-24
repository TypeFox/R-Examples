.aidsCoefNamesEst <- function( nGoods, nShifter, hom, sym ) {
   result <- NULL
   for( i in 1:( nGoods - 1 ) ) {
      result <- c( result, paste( "alpha", i ),
         paste( "beta", i ) )
      start <- ifelse( sym , i, 1 )
      stop  <- ifelse( hom , nGoods - 1, nGoods )
      for( j in start:stop ){
         result <- c( result, paste( "gamma", i, j ) )
      }
      if( nShifter > 0 ) {
         for( j in 1:nShifter ){
            result <- c( result, paste( "delta", i, j ) )
         }
      }
   }
   return( result )
}

.aidsCoefNamesAll <- function( nGoods, nShifter ) {
   result <- c(
      paste( "alpha", c( 1:nGoods ) ),
      paste( "beta", c( 1:nGoods ) ),
      paste( "gamma", rep( 1:nGoods, each = nGoods ),
         rep( 1:nGoods, nGoods ) ) )
   if( nShifter > 0 ) {
      result <- c( result,
         paste( "delta", rep( 1:nGoods, each = nShifter ),
         rep( 1:nShifter, nGoods ) ) )
   }
   return( result )
}
