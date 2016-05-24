## ===== calculation of elasticities from beta matrix ===
snqProfitElaJacobian <- function( beta, prices, quant, weights ) {
   if( !is.matrix( beta ) ) {
      stop( "argument 'beta' must be a matrix" )
   }
   if( nrow( beta ) != ncol( beta ) ) {
      stop( "argument 'beta' must be a quadratic matrix" )
   }
   if( length( prices ) != length( quant ) ) {
      stop( "arguments 'prices' and 'quant' must have the same length" )
   }
   if( length( prices ) != length( weights ) ) {
      stop( "arguments 'prices' and 'weights' must have the same length" )
   }
   if( nrow( beta ) != length( prices ) ) {
      stop( "arguments 'prices' must have as many elements as",
         " argument 'beta' has rows" )
   }
   nNetput  <- ncol( beta )
   prices   <- unlist( prices )
   quant    <- unlist( quant )
   normPrice <- sum( t( prices ) %*% weights )
   quantNames   <- .snqProfitQuantNames( quant, nNetput )
   priceNames   <- .snqProfitPriceNames( prices, nNetput )

   jacobian <- matrix( 0, nrow = nNetput^2, ncol = nNetput^2 )
   rownames( jacobian ) <- paste( "E", rep( quantNames, each = nNetput ),
      rep( priceNames, nNetput ) )
   colnames( jacobian ) <- paste( "beta", rep( 1:nNetput, each = nNetput ),
      rep( 1:nNetput, nNetput ) )
   bName <- array( paste( "beta", rep( 1:nNetput, nNetput ),
      rep( 1:nNetput, each = nNetput ) ), dim = c( nNetput, nNetput ) )
   eName <- array( paste( "E", rep( quantNames, nNetput ),
      rep( priceNames, each = nNetput ) ), dim = c( nNetput, nNetput ) )


   for( i in 1:nNetput ) {
      for( j in 1:nNetput ) {
         jacobian[ eName[ i, j ], bName[ i, j ] ] <-
            prices[ j ] / ( quant[ i ] * normPrice )
         for( k in 1:nNetput ) {
            jacobian[ eName[ i, j ], bName[ i, k ] ] <-
               jacobian[ eName[ i, j ], bName[ i, k ] ] -
               weights[ j ] * prices[ k ] * prices[ j ] /
               ( quant[ i ] * normPrice^2 )
         }
         for( l in 1:nNetput ) {
            jacobian[ eName[ i, j ], bName[ j, l ] ] <-
               jacobian[ eName[ i, j ], bName[ j, l ] ] -
               weights[ i ] * prices[ l ] * prices[ j ] /
               ( quant[ i ] * normPrice^2 )
         }
         for( k in 1:nNetput ) {
            for( l in 1:nNetput ) {
               jacobian[ eName[ i, j ], bName[ k, l ] ] <-
                  jacobian[ eName[ i, j ], bName[ k, l ] ] +
                  weights[ i ] * weights[ j ] * prices[ k ] * prices[ l ] *
                  prices[ j ] / ( quant[ i ] * normPrice^3 )
            }
         }
      }
   }
#    # test: compare with  results of snqProfitHessianDeriv(  )
#    hessianDeriv <- jacobian
#    for( i in 1:nNetput ) {
#       for( j in 1:nNetput ) {
#          hessianDeriv[ eName[ i, j ], ] <- hessianDeriv[ eName[ i, j ], ] *
#             quant[ i ] / prices[ j ]
#       }
#    }
#    hessianDeriv[ , "beta 1 1" ] <- hessianDeriv[ , "beta 1 1" ] -
#       hessianDeriv[ , "beta 1 3" ] - hessianDeriv[ , "beta 3 1" ] +
#       hessianDeriv[ , "beta 3 3" ]
#    hessianDeriv[ , "beta 1 2" ] <- hessianDeriv[ , "beta 1 2" ] -
#       hessianDeriv[ , "beta 1 3" ] + hessianDeriv[ , "beta 2 1" ] -
#       hessianDeriv[ , "beta 2 3" ] - hessianDeriv[ , "beta 3 1" ] -
#       hessianDeriv[ , "beta 3 2" ] + 2 * hessianDeriv[ , "beta 3 3" ]
#    hessianDeriv[ , "beta 2 2" ] <- hessianDeriv[ , "beta 2 2" ] -
#       hessianDeriv[ , "beta 2 3" ] - hessianDeriv[ , "beta 3 2" ] +
#       hessianDeriv[ , "beta 3 3" ]
#    return( hessianDeriv )
   return( jacobian )
}
