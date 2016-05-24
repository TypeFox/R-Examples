## ---- snqProfit: default weights for normalizing prices ----------
snqProfitWeights <- function( priceNames, quantNames, data, method = "DW92", base = 1 ) {

   if( length( quantNames ) != length( priceNames ) ) {
      stop( "arguments 'quantNames' and 'priceNames' must have the same length" )
   }

   nNetput <- length( quantNames )
   weights <- array( 0, c( nNetput ) )   # weights of the netput prices

   if( method == "DW92" ) {
      totalValue <- 0
      for( i in 1:nNetput ) {
         totalValue <- totalValue + mean( abs( data[[ quantNames[ i ] ]] ) ) *
            mean( data[[ priceNames[ i ] ]][ base ] )
      }
      for( i in 1:nNetput ) {
         weights[ i ] <- mean( abs( data[[ quantNames[ i ] ]] ) ) *
            mean( data[[ priceNames[ i ] ]][ base ] ) / totalValue
      }
   } else if( method == "Kohli" ) {
      totalValue <- 0
      for( i in 1:nNetput ) {
         totalValue <- totalValue + mean( data[[ priceNames[ i ] ]] ) *
            mean( abs( data[[ quantNames[ i ] ]] ) )
      }
      for( i in 1:nNetput ) {
         weights[ i ] <-
            mean( abs( data[[ quantNames[ i ] ]] ) ) /
            totalValue
      }
   } else if( method == "Ooms" ) {
      totalValues <- 0
      for( i in 1:nNetput ) {
         totalValues <- totalValues +
            abs( data[[ priceNames[ i ] ]] * data[[ quantNames[ i ] ]] )
      }
      for( i in 1:nNetput ) {
         weights[ i ] <- mean( abs( data[[ priceNames[ i ] ]] *
            data[[ quantNames[ i ] ]] ) ) / mean( totalValues )
      }
   } else {
      stop( "argument 'method' must be either 'DW92', 'Kohli' or 'Ooms'" )
   }

   return( weights )
}
