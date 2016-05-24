.aidsJacobian <- function( allCoef, priceNames, totExpName, data = NULL,
      shifterNames = NULL, omitLast = TRUE, alpha0 = 0 ) {
   nObs <- nrow( data )
   nGoods <- length( priceNames )
   nShifter <- length( shifterNames )
   nExogEq <- 2 + nGoods + nShifter
   coef <- .aidsCoef( allCoef, nGoods = nGoods, nShifter = nShifter, alpha0 = alpha0 )
   hom <- all.equal( rowSums( coef$gamma ), rep( 0, nGoods ) ) == TRUE
   sym <- all.equal( coef$gamma, t( coef$gamma ) ) == TRUE
   lnp <- aidsPx( "TL", priceNames, coef = coef, data = data )
   result <- matrix( 0, nrow = nObs * ( nGoods - 1 ),
      ncol = nExogEq * ( nGoods - 1 ) )
   colnames( result ) <- .aidsCoefNamesEst( nGoods, nShifter,
      hom = FALSE, sym = FALSE )
   aName <- paste( "alpha", c( 1:nGoods ) )
   bName <- paste( "beta", c( 1:nGoods ) )
   gName <- matrix( paste( "gamma", rep( 1:nGoods, nGoods ),
      rep( 1:nGoods, each = nGoods ) ), nrow = nGoods, ncol = nGoods )
   if( nShifter > 0 ) {
      dName <- matrix( paste( "delta", rep( 1:nGoods, nShifter ),
         rep( 1:nShifter, each = nGoods ) ), nrow = nGoods, ncol = nShifter )
   }

   for( eq in 1:( nGoods - 1 ) ) {
      myRows <- ( ( eq - 1 ) * nObs + 1 ):( eq * nObs )
      # derivatives of alphas
      for( i in 1:( nGoods - 1 ) ) {
         result[ myRows, aName[ i ] ] <- ( i == eq ) -
            coef$beta[ eq ] *
            ( log( data[[ priceNames[ i ] ]] ) -
            log( data[[ priceNames[ nGoods ] ]] ) )
      }
      # derivatives of betas
      result[ myRows, bName[ eq ] ] <- log( data[[ totExpName ]] ) - lnp
      # derivatives of gammas
      for( i in 1:( nGoods - 1 ) ) {
         for( j in 1:( nGoods - hom ) ) {
            result[ myRows, gName[ i, j ] ] <-
               ( i == eq ) * ( log( data[[ priceNames[ j ] ]] ) -
                  hom * log( data[[ priceNames[ nGoods ] ]] ) ) -
               0.5 * coef$beta[ eq ] *
               ( log( data[[ priceNames[ i ] ]] ) -
                  log( data[[ priceNames[ nGoods ] ]] ) ) *
               ( log( data[[ priceNames[ j ] ]] ) -
                  hom * log( data[[ priceNames[ nGoods ] ]] ) )
         }
      }
      # derivatives of deltas
      if( nShifter > 0 ) {
         for( i in 1:nShifter ) {
            result[ myRows, dName[ eq, i ] ] <-
               data[[ shifterNames[ i ] ]]
         }
      }
   }
   delCols <- NULL
   for( i in 1:( nGoods - 1 ) ) {
      if( hom ) {
         delCols <- c( delCols, gName[ i, nGoods ] )
      }
      if( sym && i >= 2 ) {
         for( j in 1:( i - 1 ) ) {
            delCol <- gName[ i, j ]
            stayCol <- gName[ j, i ]
            result[ , stayCol ] <- result[ , stayCol ] + result[ , delCol ]
            delCols <- c( delCols, delCol )
         }
      }
   }
   if( !is.null( delCols ) ) {
      result <- result[ , !colnames( result ) %in% delCols ]
   }
   return( result )
}
