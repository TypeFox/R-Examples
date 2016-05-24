aidsElas <- function( coef, prices = NULL, shares = NULL, totExp = NULL,
   method = "AIDS", priceIndex = "TL", basePrices = NULL, baseShares = NULL,
   quantNames = NULL, priceNames = NULL, coefCov = NULL, df = NULL ) {

   if( !is.null( coef$delta ) ) {
      stop( "calculating demand elasticities for models with demand shifters",
         " has not been implemented yet" )
   }

   nGoods <- length( coef$alpha )

   coefCheckResult <- .aidsCheckCoef( coef, variables = list(
      list( ifelse( is.null( prices ), NA, length( prices ) ), "prices", "goods" ),
      list( ifelse( is.null( shares ), NA, length( shares ) ), "shares", "goods" ),
      list( ifelse( is.null( quantNames ), NA, length( quantNames ) ), 
         "quantNames", "goods" ),
      list( ifelse( is.null( priceNames ), NA, length( priceNames ) ), 
         "priceNames", "goods" ) ) )
   if( !is.null( coefCheckResult ) ){
      stop( coefCheckResult )
   }

   if( priceIndex == "TL" && method != "AIDS" ) {
      stop( "there is no formula/method '", method,
         "' for calculating elasticities",
         " the original AIDS (with the translog price index 'TL')" )
   }

   if( is.null( shares ) ) {
      if( is.null( prices ) ) {
         stop( "either argument 'prices' or argument 'shares'",
            " must be specified" )
      }
      if( is.null( totExp ) ) {
         stop( "if argument 'shares' is not specified,",
            " argument 'totExp' must be specified" )
      }
      tempData <- data.frame( totExp = totExp )
      tempPriceNames <- paste ( "p", c( 1:nGoods ), sep = "" )
      for( i in  1:nGoods ) {
         tempData[[ tempPriceNames[ i ] ]] <- prices[ i ]
      }
      if( priceIndex == "SL" ) {
         tempPriceIndex <- "S"
      } else {
         tempPriceIndex <- priceIndex
      }
      shares <- as.numeric( aidsCalc( priceNames = tempPriceNames,
         totExpName = "totExp", coef = coef, data = tempData,
         priceIndex = tempPriceIndex, basePrices = basePrices,
         baseShares = baseShares )$shares )
      rm( tempData, tempPriceNames, tempPriceIndex )
   }

   if( is.null( quantNames ) ) {
      quantNames <- .aidsQuantNames( shares, coef, nGoods )
   }
   if( is.null( priceNames ) ) {
      priceNames <- .aidsPriceNames( prices, coef, nGoods )
   }

   if( method %in% c( "AIDS", "GA", "B1", "B2" ) ) {
      if( is.null( prices ) ) {
         stop( "methods 'AIDS', 'GA', 'B1', and 'B2'",
            " require argument 'prices'" )
      }
   }

   ela <- list()
   ela$method <- method
   ela$priceIndex <- priceIndex

   ones <- rep( 1, nGoods )

   if( method == "AIDS" ) {
      ela$exp <- ones + coef$beta/shares
      ela$marshall <- -diag( 1, nGoods, nGoods ) +
         coef$gamma / ( shares %*% t( ones ) ) -
         coef$beta %*% t( ones ) *
         ( ones %*% t( coef$alpha ) +
         ones %*% t( coef$gamma %*% log( prices ))) /
         ( shares %*% t( ones ) )
   } else if( method %in% c( "Ch", "Go" ) &&
         priceIndex %in% c( "S", "SL", "P" ) ) {
      ela$exp <- ones + coef$beta / shares
      ela$marshall <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( shares %*% t( ones ) ) -
         coef$beta %*% t( ones ) *
         ones %*% t( shares ) /
         ( shares %*% t( ones ) )
   } else if( method %in% c( "Go", "Ch", "B1", "GA", "B2" ) &&
         priceIndex %in% c( "L", "Ls" ) ) {
      if( is.null( baseShares ) ) {
         stop( "calculation of demand elasticities with method",
            " 'B1'/'GA' or 'B2'",
            " for models with Laspeyres or simplified Laspeyres price index",
            " requires argument 'baseShares'" )
      }
      ela$exp <- ones + coef$beta / shares
      ela$marshall <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( shares %*% t( ones ) ) -
         coef$beta %*% t( ones ) *
         ones %*% t( baseShares ) /
         ( shares %*% t( ones ) )
   } else if( method %in% c( "Ch", "Go" ) && priceIndex == "T" ) {
      ela$exp <- ones + coef$beta / shares
      ela$marshall <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( shares %*% t( ones ) ) -
         coef$beta %*% t( ones ) *
         ( ones %*% t( shares + baseShares ) ) /
         ( 2 * shares %*% t( ones ) )
   } else if( method == "EU" &&
         priceIndex %in% c( "S", "SL", "P", "L", "Ls", "T" ) ) {
      ela$exp <- ones + coef$beta / shares
      ela$marshall <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( shares %*% t( ones ) )
   } else if( method %in% c( "GA", "B1" ) && priceIndex %in% c( "S", "SL" ) ) {
      denom <- 1  # part of denominator for exp. + Marsh. elasticities
      for( k in 1:nGoods ) {
         denom <- denom + coef$beta[ k ] * log( prices[ k ] )
      }
      ela$exp <- ones + coef$beta / ( shares * denom )
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            numer <- shares[ j ] # part of numerator for Marsh. elasticities
            for( k in 1:nGoods ) {
               numer <- numer + coef$gamma[ k, j ] * log( prices[ k ] )
            }
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / shares[ i ] -
               ( coef$beta[ i ] / shares[ i ] ) * numer / denom
         }
      }
   } else if( method %in% c( "GA", "B1" ) && priceIndex %in% c( "P" ) ) {
      if( is.null( basePrices ) ) {
         stop( "calculations of demand elasticities with method 'B1'/'GA'",
            " for models with Paasche price index",
            " require argument 'basePrices'" )
      }
      denom <- 1  # part of denominator for exp. + Marsh. elasticities
      for( k in 1:nGoods ) {
         denom <- denom + coef$beta[ k ] * log( prices[ k ] / basePrices[ k ] )
      }
      ela$exp <- ones + coef$beta / ( shares * denom )
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            numer <- shares[ j ] # part of numerator for Marsh. elasticities
            for( k in 1:nGoods ) {
               numer <- numer + coef$gamma[ k, j ] *
                  log( prices[ k ] / basePrices[ k ] )
            }
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / shares[ i ] -
               ( coef$beta[ i ] / shares[ i ] ) * numer / denom
         }
      }
   } else if( method %in% c( "GA", "B1" ) && priceIndex %in% c( "T" ) ) {
      if( is.null( basePrices ) ) {
         stop( "calculations of demand elasticities with method 'B1'/'GA'",
            " for models with Tornqvist price index",
            " require argument 'basePrices'" )
      }
      if( is.null( baseShares ) ) {
         stop( "calculations of demand elasticities with method 'B1'/'GA'",
            " for models with Tornqvist price index",
            " require argument 'baseShares'" )
      }
      denom <- 1  # part of denominator for exp. + Marsh. elasticities
      for( k in 1:nGoods ) {
         denom <- denom + 0.5 * coef$beta[ k ] *
            log( prices[ k ] / basePrices[ k ] )
      }
      ela$exp <- ones + coef$beta / ( shares * denom )
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            numer <- shares[ j ] + baseShares[ j ]
               # part of numerator for Marsh. elasticities
            for( k in 1:nGoods ) {
               numer <- numer + coef$gamma[ k, j ] *
                  log( prices[ k ] / basePrices[ k ] )
            }
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / shares[ i ] -
               ( coef$beta[ i ] / ( 2 * shares[ i ] ) ) * numer / denom
         }
      }
   } else if( method %in% c( "B2" ) && priceIndex %in% c( "S", "SL" ) ) {
      paren <- 1  # term in parenthesis for expenditure elasticities
      for( k in 1:nGoods ) {
         paren <- paren - coef$beta[ k ] * log( prices[ k ] )
      }
      ela$exp <- ones + ( coef$beta / shares ) * paren
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            parenSmall <- coef$alpha[ j ] # term in small par. for Marsh. elast.
            for( l in 1:nGoods ) {
               parenSmall <- parenSmall +
                  coef$gamma[ l, j ] * log( prices[ l ] )
            }
            parenBig <- shares[ j ] # term in big parenthesis for Marsh. elast.
            for( k in 1:nGoods ) {
               parenBig <- parenBig +
                  coef$gamma[ k, j ] * log( prices[ k ] ) -
                  coef$beta[ k ] * log( prices[ k ] ) * parenSmall
            }
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / shares[ i ] -
               ( coef$beta[ i ] / shares[ i ] ) * parenBig
         }
      }
   } else if( method %in% c( "B2" ) && priceIndex == "P" ) {
      if( is.null( basePrices ) ) {
         stop( "calculations of demand elasticities with method 'B2'",
            " for models with Paasche price index",
            " require argument 'basePrices'" )
      }
      paren <- 1  # term in parenthesis for expenditure elasticities
      for( k in 1:nGoods ) {
         paren <- paren - coef$beta[ k ] * log( prices[ k ] / basePrices[ k ] )
      }
      ela$exp <- ones + ( coef$beta / shares ) * paren
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            parenSmall <- coef$alpha[ j ] # term in small par. for Marsh. elast.
            for( l in 1:nGoods ) {
               parenSmall <- parenSmall +
                  coef$gamma[ l, j ] * log( prices[ l ] )
            }
            parenBig <- shares[ j ] # term in big parenthesis for Marsh. elast.
            for( k in 1:nGoods ) {
               parenBig <- parenBig +
                  coef$gamma[ k, j ] * log( prices[ k ] / basePrices[ k ] ) -
                  coef$beta[ k ] * log( prices[ k ] / basePrices[ k ] ) *
                  parenSmall
            }
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / shares[ i ] -
               ( coef$beta[ i ] / shares[ i ] ) * parenBig
         }
      }
   } else if( method %in% c( "B2" ) && priceIndex == "T" ) {
      if( is.null( basePrices ) ) {
         stop( "calculations of demand elasticities with method 'B2'",
            " for models with Tornqvist price index",
            " require argument 'basePrices'" )
      }
      if( is.null( baseShares ) ) {
         stop( "calculations of demand elasticities with method 'B2'",
            " for models with Tornqvist price index",
            " require argument 'baseShares'" )
      }
      paren <- 1  # term in parenthesis for expenditure elasticities
      for( k in 1:nGoods ) {
         paren <- paren - 0.5 * coef$beta[ k ] *
            log( prices[ k ] / basePrices[ k ] )
      }
      ela$exp <- ones + ( coef$beta / shares ) * paren
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            parenSmall <- coef$alpha[ j ] # term in small par. for Marsh. elast.
            for( l in 1:nGoods ) {
               parenSmall <- parenSmall +
                  coef$gamma[ l, j ] * log( prices[ l ] )
            }
            parenBig <- shares[ j ] + baseShares[ j ]
               # term in big parenthesis for Marsh. elast.
            for( k in 1:nGoods ) {
               parenBig <- parenBig +
                  coef$gamma[ k, j ] * log( prices[ k ] / basePrices[ k ] ) -
                  coef$beta[ k ] * log( prices[ k ] / basePrices[ k ] ) *
                  parenSmall
            }
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / shares[ i ] -
               ( coef$beta[ i ] / ( 2 * shares[ i ] ) ) * parenBig
         }
      }
   } else {
      stop( "calculations of demand elasticities",
         " with elasticity formula '", method, "' (argument 'method')",
         " for models",
         " with price index '", priceIndex, "' (argument 'priceIndex')",
         " are currently not (yet) implemented" )
   }
   ela$hicks <- ela$marshall + ( ela$exp %*% t( ones ) ) *
      ( ones %*% t(shares))
   names( ela$exp )         <- quantNames
   rownames( ela$hicks )    <- quantNames
   colnames( ela$hicks )    <- priceNames
   rownames( ela$marshall ) <- quantNames
   colnames( ela$marshall ) <- priceNames
   if( !is.null( coefCov ) && method %in% c( "AIDS" ) ) {
      jacobian <- .aidsElasJacobian( coef = coef, shares = shares, prices = prices,
         method = method, quantNames = quantNames, priceNames = priceNames )
      ela$allVcov      <- jacobian$all      %*% coefCov %*% t( jacobian$all )
      ela$expVcov      <- jacobian$exp      %*% coefCov %*% t( jacobian$exp )
      ela$hicksVcov    <- jacobian$hicks    %*% coefCov %*% t( jacobian$hicks )
      ela$marshallVcov <- jacobian$marshall %*% coefCov %*% t( jacobian$marshall )
      # standard errors
      ela$expStEr      <- diag( ela$expVcov )^0.5
      ela$hicksStEr    <- matrix( diag( ela$hicksVcov )^0.5,
         ncol = nGoods, byrow = TRUE )
      ela$marshallStEr <-  matrix( diag( ela$marshallVcov )^0.5,
         ncol = nGoods, byrow = TRUE )
      # dim names for standard errors
      names( ela$expStEr )         <- names( ela$exp )
      dimnames( ela$hicksStEr )    <- dimnames( ela$hicks )
      dimnames( ela$marshallStEr ) <- dimnames( ela$marshall )
      # t-values
      ela$expTval      <- ela$exp      / ela$expStEr
      ela$hicksTval    <- ela$hicks    / ela$hicksStEr
      ela$marshallTval <- ela$marshall / ela$marshallStEr
      if( !is.null( df ) ) {
         ela$expPval <- 2 * pt( abs( ela$expTval ), df,
            lower.tail = FALSE )
         ela$hicksPval <- 2 * pt( abs( ela$hicksTval ), df,
            lower.tail = FALSE )
         ela$marshallPval <- 2 * pt( abs( ela$marshallTval ), df,
            lower.tail = FALSE )
      }
      ela$df <- df
   }
   class( ela ) <- "aidsElas"
   return( ela )
}
