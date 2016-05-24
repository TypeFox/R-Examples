.aidsCoef <- function( coef, nGoods, nShifter = 0, cov = NULL, df = 1,
      LA = TRUE, priceNames = NULL, shareNames = NULL, shifterNames = NULL,
      alpha0 = NULL ) {
   # nGoods <- -0.5 + ( 2.25 + nrow( array( coef ) ) )^0.5
   nExogEq <- nGoods + 2 + nShifter
   if( LA ) {
      M <- matrix( 0, nGoods * nExogEq, ( nGoods - 1 ) * nExogEq )
      rownames( M ) <- .aidsCoefNamesAll( nGoods, nShifter )
      colnames( M ) <- .aidsCoefNamesEst( nGoods, nShifter, hom = FALSE,
         sym = FALSE )
      aName <- paste( "alpha", c( 1:nGoods ) )
      bName <- paste( "beta", c( 1:nGoods ) )
      gName <- matrix( paste( "gamma", rep( 1:nGoods, nGoods ),
         rep( 1:nGoods, each = nGoods ) ), nrow = nGoods, ncol = nGoods )
      if( nShifter > 0 ) {
         dName <- matrix( paste( "delta", rep( 1:nGoods, nShifter ),
            rep( 1:nShifter, each = nGoods ) ), nrow = nGoods, ncol = nShifter )
      }
      for(i in 1:( nGoods - 1 ) ) {
         M[ aName[ i ], aName[ i ] ]   <-  1   # alpha[i]
         M[ aName[ nGoods ], aName[ i ] ]   <- -1   # alpha[nGoods]
         M[ bName[ i ], bName[ i ] ] <-  1   # beta[i]
         M[ bName[ nGoods ], bName[ i ] ] <- -1 # beta[ nGoods ]
         for( j in 1:nGoods ) {
            M[ gName[ i, j ], gName[ i, j ] ] <-  1   # gamma[i,j]
            M[ gName[ nGoods, j ], gName[ i, j ] ] <- -1   # gamma[nGoods,j]
         }
        if( nShifter > 0 ){
            for( j in 1:nShifter ) {
               M[ dName[ i, j ], dName[ i, j ] ] <-  1   # gamma[i,j]
               M[ dName[ nGoods, j ], dName[ i, j ] ] <- -1   # gamma[nGoods,j]
            }
         }
      }
      all     <- c( M %*% coef )
      names( all ) <- .aidsCoefNamesAll( nGoods, nShifter )
      all[ aName[ nGoods ] ]  <- all[ aName[ nGoods ] ] + 1
      alpha   <- all[1:nGoods]
      beta    <- all[(nGoods+1):(2*nGoods)]
      gamma   <- matrix( all[(2*nGoods+1):(nGoods*(nGoods+2))],
         nrow = nGoods, ncol = nGoods, byrow = TRUE )
      if( nShifter > 0 ){
         delta <- matrix( all[ ( nGoods * ( nGoods + 2 ) + 1 ):length( all ) ],
            nrow = nGoods, ncol = nShifter, byrow = TRUE )
      } else {
         delta <- NULL
      }
      allcov <- NULL
      stat    <- NULL
      if(!is.null(cov)) {
         allcov   <- M %*% cov %*% t(M)
         rownames( allcov ) <- names( all )
         colnames( allcov ) <- names( all )
         stat     <- coefTable( all, sqrt( diag( allcov ) ), df )
      }
   }
   if( !is.null( shareNames ) ) {
      names( alpha ) <- shareNames
      names( beta ) <- shareNames
      rownames( gamma ) <- shareNames
      if( !is.null( delta ) ){
         rownames( delta ) <- shareNames
      }
   }
   if( !is.null( priceNames ) ) {
      colnames( gamma ) <- priceNames
   }
   if( !is.null( delta ) & !is.null( shifterNames ) ){
      colnames( delta ) <- shifterNames
   }
   result <- list()
   result$alpha0 <- alpha0
   result$alpha <- alpha
   result$beta  <- beta
   result$gamma <- gamma
   result$delta <- delta
   result$all   <- all
   result$allcov <- allcov
   result$stat  <- stat
   return( result )
}
