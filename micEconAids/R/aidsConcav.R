aidsConcav <- function( priceNames, totExpName, coef, data,
      shareNames = NULL ) {

   if( !is.null( shareNames ) && length( priceNames ) != length( shareNames ) ) {
      stop( "arguments 'priceNames' and 'shareNames' must have the same length" )
   }
   if( is.null( coef$alpha0 ) ) {
      stop( "argument 'coef' must have element 'alpha0'" )
   }

   if( isSymmetric( coef$gamma, tol = 1e-6, check.attributes = FALSE )
         != TRUE ) {
      stop( "there does not exist an expenditure function,",
         " because the matrix of 'gamma' coefficients is not symmetric" )
   }

   result <- list()
   nGoods <- length( priceNames )
   nObs <- nrow( data )

   xt <- data[[ totExpName ]]
   shareMat <- array( NA, c( nObs, nGoods ) )
   for( i in 1: nGoods ) {
      if( !is.null( shareNames ) ) {
         shareMat[ , i ] <- data[[ shareNames[ i ] ]]
      }
   }
   fitted <- aidsCalc( priceNames, totExpName, data = data,
      coef = coef )
   if( is.null( shareNames ) ) {
      shareMat <- as.matrix( fitted$shares )
   }

   # checking concavity
   result$cMatrices <- list()
   result$concavity <- rep( NA, nObs )

   lnp <- aidsPx( "TL", priceNames, data = data, coef = coef )

   for( t in 1:nObs ) {
      result$cMatrices[[ t ]] <- coef$gamma + ( coef$beta %*% t( coef$beta ) ) *
         ( log( xt[ t ] ) - lnp[ t ] ) -
         diag( shareMat[ t, ] ) + shareMat[ t, ] %*% t( shareMat[ t, ] )

      result$concavity[ t ] <- semidefiniteness( result$cMatrices[[ t ]][ 1:( nGoods - 1),
         1:( nGoods - 1) ], positive = FALSE )
   }

   result$nValidObs <- sum( !is.na( result$concavity ) )
   result$nConcavObs <- sum( result$concavity, na.rm = TRUE )
   result$concavPercent <- 100 * result$nConcavObs / result$nValidObs

   class( result ) <- "aidsConcav"
   return( result )
}
