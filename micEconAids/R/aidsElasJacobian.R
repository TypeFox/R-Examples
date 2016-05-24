.aidsElasJacobian <- function( coef, shares, prices = NULL, method = "AIDS",
      quantNames = NULL, priceNames = NULL ) {

   nGoods <- length( coef$alpha )
   nCoef  <- ( nGoods + 2 ) * nGoods

   if( length( coef$alpha ) != length( coef$beta ) ) {
      stop( "arguments 'alpha' and 'beta' must have the same length" )
   } else if( nrow( coef$gamma ) != ncol( coef$gamma ) ) {
      stop( "argument 'gamma' must be a square matrix" )
   } else if( length( coef$alpha ) != nrow( coef$gamma ) ) {
      stop( "number of rows of argument 'gamma' must be equal",
         " to the length of argument 'alpha'" )
   } else if(  length( coef$alpha ) != length( shares ) ) {
      stop( "arguments 'alpha' and 'shares' must have the same length" )
   } else if(  length( coef$alpha ) != length( prices ) && !is.null( prices ) ) {
      stop( "arguments 'alpha' and 'prices' must have the same length" )
   }
   if( is.null( quantNames ) ) {
      quantNames <- .aidsQuantNames( shares, coef, nGoods )
   } else {
      if( length( quantNames ) != nGoods ) {
         stop( "argument 'quantNames' must have ", nGoods, " elements" )
      }
   }
   if( is.null( priceNames ) ) {
      priceNames <- .aidsPriceNames( prices, coef, nGoods )
   } else {
      if( length( priceNames ) != nGoods ) {
         stop( "argument 'priceNames' must have ", nGoods, " elements" )
      }
   }

   if( method %in% c( "AIDS" ) ) {
      if( is.null( prices ) ) {
         stop( "the 'AIDS' method requires argument 'prices'" )
      }
   }

   createMatrix <- function( nGoods, nCoef, dim, symbol, quantNames, priceNames ) {
      result <- matrix( 0, nrow = nGoods^dim, ncol = nCoef )
      if( dim == 1 ) {
         rownames( result ) <- paste( symbol, quantNames )
      } else if( dim == 2 ) {
         rownames( result ) <- paste( symbol, rep( quantNames, each = nGoods ),
            rep( priceNames, nGoods ) )
      }
      colnames( result ) <- .aidsCoefNamesAll( nGoods, 0 )
      return( result )
   }

   jacobian <- list()
   jacobian$method   <- method
   jacobian$exp      <- createMatrix( nGoods, nCoef, 1, "Ex", quantNames, priceNames )
   jacobian$hicks    <- createMatrix( nGoods, nCoef, 2, "Eh", quantNames, priceNames )
   jacobian$marshall <- createMatrix( nGoods, nCoef, 2, "Em", quantNames, priceNames )

   shares <- array( shares )

   aName <- paste( "alpha", c( 1:nGoods ) )
   bName <- paste( "beta", c( 1:nGoods ) )
   gName <- array( paste( "gamma", rep( 1:nGoods, nGoods ),
      rep( 1:nGoods, each = nGoods ) ), dim = c( nGoods, nGoods ) )
   ehName <- array( paste( "Eh", rep( quantNames, nGoods ),
      rep( priceNames, each = nGoods ) ), dim = c( nGoods, nGoods ) )
   emName <- array( paste( "Em", rep( quantNames, nGoods ),
      rep( priceNames, each = nGoods ) ), dim = c( nGoods, nGoods ) )

   if( method == "AIDS" ) {
      prices <- array( prices )
      for( i in 1:nGoods ) {
         # expenditure elasticities
         jacobian$exp[ paste( "Ex", quantNames[ i ] ), bName[ i ] ] <-
            1 / shares[ i ]
         for( j in 1:nGoods ) {
            # Hicksian price elasticities
            jacobian$hicks[ ehName[ i, j ], aName[ j ] ] <-
               -coef$beta[ i ] / shares[ i ]
            jacobian$hicks[ ehName[ i, j ], bName[ i ] ] <-
               - ( coef$alpha[ j ] - shares[ j ] +
               coef$gamma[ j , ] %*% log( prices ) ) / shares[ i ]
            for( k in 1:nGoods ) {
               jacobian$hicks[ ehName[ i, j ], gName[ k, j ] ] <-
                  ( i == k ) / shares[ i ] -
                  coef$beta[ i ] * log( prices[ k ] ) / shares[ i ]
            }
            # Marshallian price elasticities
            jacobian$marshall[ emName[ i, j ], aName[ j ] ] <-
               -coef$beta[ i ] / shares[ i ]
            jacobian$marshall[ emName[ i, j ], bName[ i ] ] <-
               - ( coef$alpha[ j ] +
               coef$gamma[ j , ] %*% log( prices ) ) / shares[ i ]
            for( k in 1:nGoods ) {
               jacobian$marshall[ emName[ i, j ], gName[ k, j ] ] <-
                  ( i == k ) / shares[ i ] -
                  coef$beta[ i ] * log( prices[ k ] ) / shares[ i ]
            }
         }
      }
   } else {
      stop( "method '", as.character( method ), "' is not supported" )
   }
   jacobian$all <- rbind( jacobian$exp, jacobian$hicks, jacobian$marshall )
   return( jacobian )
}
