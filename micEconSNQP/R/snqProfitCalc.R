snqProfitCalc <- function( priceNames, fixNames, data, weights,
      scalingFactors = rep( 1, length( weights ) ), coef,
      quantNames = NULL, form = 0 ) {

   checkNames( c( priceNames, fixNames ), names( data ) )

   nNetput <- length( priceNames )
   nFix    <- length( fixNames )
   nObs    <- nrow( data )

   snqProfitTestCoef( nNetput, nFix, coef, form = form )

   if( is.null( quantNames ) ) {
      quantNames <- paste( "X", 1:nNetput, sep = "" )
   }
   normPrice <- numeric( nObs )
   for( i in 1:nNetput ) {
      normPrice <- normPrice + data[[ priceNames[ i ] ]] * scalingFactors[ i ] *
         weights[ i ]
   }

   qNetput <- array( 0, c( nObs, nNetput ) )
   for( i in 1:nNetput ) {
      qNetput[ , i ] <- coef$alpha[ i ]
      for( j in 1:nNetput ) {
         qNetput[ , i ] <- qNetput[ , i ] +
            coef$beta[ i, j ] * data[[ priceNames[ j ] ]] * scalingFactors[ j ] /
            normPrice
         for( k in 1:nNetput ) {
            qNetput[ , i ] <- qNetput[ , i ] -
               0.5 * weights[ i ] * coef$beta[ j, k ] *
               data[[ priceNames[ j ] ]] * scalingFactors[ j ] *
               data[[ priceNames[ k ] ]] * scalingFactors[ k ] / normPrice^2
         }
      }
      if( nFix > 0 ) {
         for( j in 1:nFix ) {
            qNetput[ , i ] <- qNetput[ , i ] +
               coef$delta[ i, j ] * data[[ fixNames[ j ] ]]
            for( k in 1:nFix ) {
               if( form == 0 ) {
                  qNetput[ , i ] <- qNetput[ , i ] + 0.5 * weights[ i ] *
                     coef$gamma[ j, k ] * data[[ fixNames[ j ] ]] * data[[ fixNames[ k ] ]]
               } else {
                  qNetput[ , i ] <- qNetput[ , i ] + 0.5 *
                     coef$gamma[ i, j, k ] * data[[ fixNames[ j ] ]] *
                     data[[ fixNames[ k ] ]]
               }
            }
         }
      }
   }
   # qNetput <- array(1,c(nObs)) %*% t(coef$alpha) + (P%*%coef$beta) /
   #    (normPrice %*% array(1,c(1,nNetput))) -
   #    0.5 * (diag(P %*% coef$beta %*% t(P)) %*% t(weights)) /
   #    ((normPrice^2) %*% array(1,c(1,nNetput))) +
   #    Z %*% t(coef$delta) +
   #    0.5 * diag(Z %*% coef$gamma %*% t(Z)) %*% t(weights)

   profit <- numeric( nObs )
   for( i in 1:nNetput ) {
      profit <- profit + coef$alpha[ i ] * data[[ priceNames[ i ] ]] *
         scalingFactors[ i ]
      for( j in 1:nNetput ) {
         profit <- profit + 0.5 * coef$beta[ i, j ] *
            data[[ priceNames[ i ] ]] * scalingFactors[ i ] *
            data[[ priceNames[ j ] ]] * scalingFactors[ j ] / normPrice
      }
   }
   if( nFix > 0 ) {
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            profit <- profit + coef$delta[ i, j ] * data[[ priceNames[ i ] ]] *
               scalingFactors[ i ] * data[[ fixNames[ j ] ]]
         }
      }
      if( form == 0 ) {
         for( i in 1:nFix ) {
            for( j in 1:nFix ) {
               profit <- profit + 0.5 * normPrice * coef$gamma[ i, j ] *
                  data[[ fixNames[ i ] ]] * data[[ fixNames[ j ] ]]
            }
         }
      } else {
         for( i in 1:nNetput ) {
            for( j in 1:nFix ) {
               for( k in 1:nFix ) {
                  profit <- profit + 0.5 * coef$gamma[ i, j, k ] *
                     data[[ priceNames[ i ] ]] * scalingFactors[ i ] *
                     data[[ fixNames[ j ] ]] * data[[ fixNames[ k ] ]]
               }
            }
         }
      }
   }
   result <- as.data.frame( cbind( qNetput, profit ) )
   names( result ) <- c( quantNames, "profit" )

   return( result )
}
