aidsCalc <- function( priceNames, totExpName, coef, data,
      priceIndex = "TL", basePrices = NULL, baseShares = NULL ) {

   # check argument 'coef' (coefficients)
   coefCheckResult <- .aidsCheckCoef( coef, variables = list(
      list( length( priceNames ), "prices", "goods"  ) ) )
   if( !is.null( coefCheckResult ) ){
      stop( coefCheckResult )
   }

   # checking argument 'data'
   if( class( data ) != "data.frame" ) {
      stop( "argument 'data' must be a data frame" )
   }

   # checking (mainly) argument 'priceIndex'
   if( is.character( priceIndex ) ) {
      if( ! priceIndex %in% c( "TL", "S", "SL", "P", "L", "Ls", "T" ) ) {
         stop( "argument 'priceIndex' must be either",
            " 'TL' (translog), 'S' (Stone), 'SL' (Stone index with lagged shares),",
            " 'P' (Paasche), 'L' (Laspeyres),",
            " 'Ls' (Laspeyres, simplified), 'T' (Tornqvist), or a numeric vector",
            " providing the log values of the price index" )
      }
      if( priceIndex == "TL" && is.null( coef$alpha0 ) ) {
         stop( "calculations with the translog (TL) price index require",
            " coefficient alpha_0 (coef$alpha0)" )
      }
   } else if( is.numeric( priceIndex ) ) {
      if( length( priceIndex ) != nrow( data ) ) {
         stop( "if argument 'priceIndex' provides the values",
            " of the log price index,",
            " it must have the same length as there are observations",
            " in argument 'data'" )
      }
   } else {
      stop( "argument 'priceIndex' must be either a character string",
         " or a numeric vector" )
   }

   # tests for arguments basePrices and baseShares
   if( is.character( priceIndex ) ) {
      # basePrices
      if( priceIndex %in% c( "P", "L", "T" ) ) {
         if( is.null( basePrices ) ) {
            stop( "calculations with Paasche ('L'), Laspeyres ('L'),",
               " or Tornqvist ('T') price index require",
               " argument 'basePrices'" )
         }
         if( ! is.numeric( basePrices ) ) {
            stop( "argument 'basePrices' must be numeric" )
         }
         if( length( basePrices ) != length( priceNames ) ) {
            stop( "arguments 'basePrices' and 'priceNames' must have",
               " the same length" )
         }
      }
      # baseShares
      if( priceIndex %in% c( "L", "Ls", "T" ) ) {
         if( is.null( baseShares ) ) {
            stop( "calculations with Laspeyres ('Ls'),",
               " simplified Laspeyres ('Ls'), or",
               " Tornqvist ('T') price index require",
               " argument 'baseShares'" )
         }
         if( ! is.numeric( baseShares ) ) {
            stop( "argument 'baseShares' must be numeric" )
         }
         if( length( baseShares ) != length( priceNames ) ) {
            stop( "arguments 'baseShares' and 'priceNames' must have",
               " the same length" )
         }
      }
   }

   if( is.character( priceIndex ) ) {
      if( priceIndex == "TL" ) {
         # calculation of translog price index
         priceIndex <- aidsPx( priceIndex, priceNames, data = data, coef = coef )
      } else if( priceIndex == "L" ) {
         # calculation of Laspeyres price index
         priceIndex <- aidsPx( priceIndex, priceNames, data = data,
            coef = coef, base = list( prices = basePrices, shares = baseShares ) )
      } else if( priceIndex == "Ls" ) {
         # calculation of simplified Laspeyres price index
         priceIndex <- aidsPx( priceIndex, priceNames, data = data,
            coef = coef, base = list( shares = baseShares ) )
      }
   }

   # number of goods
   nGoods <- length( priceNames )

   shareData <- as.data.frame( matrix( NA, nrow = nrow( data ), ncol = nGoods ) )
   names( shareData ) <- paste( "w", as.character( 1:nGoods ), sep = "" )
   rownames( shareData ) <- rownames( data )
   quant <- as.data.frame( matrix( 0, nrow = nrow( data ), ncol = nGoods ) )
   names( quant ) <- paste( "q", as.character( 1:nGoods ), sep = "" )
   rownames( quant ) <- rownames( data )
   if( is.numeric( priceIndex ) ) {
      for( i in 1:nGoods ) {
         shareData[ , i ] <- coef$alpha[ i ] + coef$beta[ i ] *
            ( log( data[[ totExpName ]] ) - priceIndex )
         for( j in 1:nGoods ) {
            shareData[ , i ] <- shareData[ , i ] + coef$gamma[ i, j ] *
               log( data[[ priceNames[ j ] ]] )
         }
      }
   } else if( priceIndex == "S" ) {
      for( i in 1:nrow( data ) ) {
         logPrices <- log( as.numeric( data[ i, priceNames ] ) )
         logTotExp <- log( data[ i, totExpName ] )
         if( all( !is.na( c( logPrices, logTotExp ) ) ) ) {
            shareData[ i, ] <-
               solve( diag( nGoods ) + coef$beta %*% t( logPrices ),
                  coef$alpha + coef$gamma %*% logPrices + coef$beta * logTotExp )
         }
      }
   } else if( priceIndex == "SL" ) {
      logPrices <- log( as.numeric( data[ 1, priceNames ] ) )
      logTotExp <- log( data[ 1, totExpName ] )
      if( all( !is.na( c( logPrices, logTotExp ) ) ) ) {
         shareData[ 1, ] <-
               solve( diag( nGoods ) + coef$beta %*% t( logPrices ),
                  coef$alpha + coef$gamma %*% logPrices + coef$beta * logTotExp )
      }
      for( i in 2:nrow( data ) ) {
         logPrices <- log( as.numeric( data[ i, priceNames ] ) )
         logTotExp <- log( data[ i, totExpName ] )
         if( all( !is.na( c( logPrices, logTotExp ) ) ) ) {
            shareData[ i, ] <-
               coef$alpha + coef$gamma %*% logPrices + coef$beta * logTotExp -
               coef$beta * t( logPrices ) %*% as.numeric( shareData[ i-1, ] )
         }
      }
   } else if( priceIndex == "P" ) {
      for( i in 1:nrow( data ) ) {
         prices <- as.numeric( data[ i, priceNames ] )
         logTotExp <- log( data[ i, totExpName ] )
         if( all( !is.na( c( prices, logTotExp ) ) ) ) {
            shareData[ i, ] <-
               solve( diag( nGoods ) + coef$beta %*% t( log( prices / basePrices ) ),
                  coef$alpha + coef$gamma %*% log( prices ) + coef$beta * logTotExp )
         }
      }
   } else if( priceIndex == "T" ) {
      for( i in 1:nrow( data ) ) {
         prices <- as.numeric( data[ i, priceNames ] )
         logTotExp <- log( data[ i, totExpName ] )
         if( all( !is.na( c( prices, logTotExp ) ) ) ) {
            shareData[ i, ] <-
               solve( diag( nGoods ) + 0.5 * coef$beta %*%
                  t( log( prices / basePrices ) ),
                  coef$alpha + coef$gamma %*% log( prices ) +
                  coef$beta * logTotExp - 0.5 * coef$beta *
                  ( t( log( prices / basePrices ) ) %*% baseShares ) )
         }
      }
   } else {
      stop( "internal error" )
   }
   for( i in 1:nGoods ) {
      quant[ , i ] <- shareData[ , i ] * data[[ totExpName ]] / data[[ priceNames[ i ] ]]
   }
   result <- list()
   result$shares <- shareData
   result$quant  <- quant
   return( result )
}
