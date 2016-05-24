snqProfitShadowPrices <- function( priceNames, fixNames, estResult = NULL,
   data = estResult$data, weights = estResult$weights,
   scalingFactors = estResult$scalingFactors,
   coef = estResult$coef, form = estResult$form ) {

   checkNames( c( priceNames, fixNames ), names( data ) )

   nNetput <- length( priceNames )
   nFix    <- length( fixNames )
   nObs    <- nrow( data )

   snqProfitTestCoef( nNetput, nFix, coef, form = form,
      coefNames = c( "delta", "gamma" ) )

   normPrice <- numeric( nObs )
   for( i in 1:nNetput ) {
      normPrice <- normPrice + data[[ priceNames[ i ] ]] * scalingFactors[ i ] *
         weights[ i ]
   }

   shadowPrices <- array( 0, c( nObs, nFix ) )
   for( j in 1:nFix ) {
      for( i in 1:nNetput ) {
         shadowPrices[ , j ] <- shadowPrices[ , j ] + coef$delta[ i, j ] *
            data[[ priceNames[ i ] ]] * scalingFactors[ i ]
      }
      if( form == 0 ) {
         for( k in 1:nFix ) {
            shadowPrices[ , j ] <- shadowPrices[ , j ] + normPrice *
               coef$gamma[ j, k ] * data[[ fixNames[ k ] ]]
         }
      } else {
         for( i in 1:nNetput ) {
            for( k in 1:nFix ) {
               shadowPrices[ , j ] <- shadowPrices[ , j ] +
                  coef$gamma[ i, j, k ] * data[[ priceNames[ i ] ]] *
                  scalingFactors[ i ] * data[[ fixNames[ k ] ]]
            }
         }
      }
   }

   result <- as.data.frame( shadowPrices )
   names( result ) <- fixNames

   return( result )
}
