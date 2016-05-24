aidsPx <- function( priceIndex, priceNames, data, shareNames = NULL, base = 1,
   coef = NULL, shifterNames = NULL ) {

   if( priceIndex == "TL" ){
      if( is.null( coef ) ) {
         stop( "argument 'coef' must be specified to calculate the translog",
            " price index" )
      } else {
         coefCheckResult <- .aidsCheckCoef( coef, variables = list(
            list( length( priceNames ), "priceNames", "goods" ),
            list( ifelse( is.null( shareNames ), NA, length( shareNames ) ), 
               "shareNames", "goods" ) ) )
         if( !is.null( coefCheckResult ) ){
            stop( coefCheckResult )
         }
         if( is.null( coef$alpha0 ) ) {
            stop( "argument 'coef' must have element 'alpha0'" )
         }
      }
   } else {
      if( is.null( shareNames ) &&
            !( priceIndex %in% c( "L", "Ls" ) && class( base ) == "list" ) ) {
         stop( "argument 'shareNames' must must be specified to calculate",
            " price index '", priceIndex, "'" )
      }
   }

   nGoods <- length( priceNames )
   nShifter <- length( shifterNames )
   nObs <- nrow(  data )
   lnp <- numeric( nObs )

   if( priceIndex %in% c( "L", "P", "T" ) ){
      if( class( base ) == "list" ){
         if( is.null( base$prices ) ){
            stop( "if argument 'priceIndex' is '", priceIndex, "'",
               " and argument 'base' is a list,",
               " this list must have an element 'prices'" )
         } else if( length( base$prices ) != nGoods ){
            stop( "element 'prices' of argument 'base'",
               " must have the same length as argument 'priceNames'" )
         } else {
            basePrices <- base$prices
         }
      } else {
         basePrices <- rep( NA, nGoods )
         for( i in 1:nGoods ) {
            basePrices[ i ] <- mean( data[[ priceNames[ i ] ]][ base ] )
         }
      }
   }
   if( priceIndex %in% c( "L", "Ls", "T" ) ){
      if( class( base ) == "list" ){
         if( is.null( base$shares ) ){
            stop( "if argument 'priceIndex' is '", priceIndex, "'",
               " and argument 'base' is a list,",
               " this list must have an element 'shares'" )
         } else if( length( base$shares ) != nGoods ){
            stop( "element 'shares' of argument 'base'",
               " must have the same length as argument 'priceNames'" )
         } else {
            if( all.equal( sum( base$shares ), 1 ) != TRUE ){
               warning( "the base expenditure shares specified",
                  " by element 'shares' of argument 'base'",
                  " do not sum up to 1 (deviation from 1 = ",
                  formatC( sum( base$shares ) - 1, digits = 3, format = "g" ),
                  ")" )
            }
            baseShares <- base$shares
         }
      } else {
         baseShares <- rep( NA, nGoods )
         for( i in 1:nGoods ) {
            baseShares[ i ] <- mean( data[[ shareNames[ i ] ]][ base ] )
         }
      }
   }

   if(priceIndex=="S") {      # Stone index
      for( i in 1:nGoods ) {
         lnp <- lnp + data[[ shareNames[ i ] ]] * log( data[[ priceNames[ i ] ]] )
      }
   } else if(priceIndex=="SL") {     # Stone index with lagged shares
      lnp[ 1 ] <- NA
      for( i in 1:nGoods ) {
         lnp[ 2:nObs ] <- lnp[ 2:nObs ] +
            data[[ shareNames[ i ] ]][ 1:(nObs-1) ] *
            log( data[[ priceNames[ i ] ]][ 2:nObs ] )
      }
   } else if(priceIndex=="P") {      # log-Paasche index
      for( i in 1:nGoods) {
         lnp <- lnp + data[[ shareNames[ i ] ]] * log( data[[ priceNames[ i ] ]] /
            basePrices[ i ] )
      }
   } else if(priceIndex=="L") {      # log-Laspeyres index
      for( i in 1:nGoods) {
         lnp <- lnp + baseShares[ i ] *
            log( data[[ priceNames[ i ] ]] / basePrices[ i ] )
      }
   } else if(priceIndex=="Ls") {      # log-Laspeyres index, simplified
      for( i in 1:nGoods) {
         lnp <- lnp + baseShares[ i ] *
            log( data[[ priceNames[ i ] ]] )
      }
   } else if(priceIndex=="T") {      # Tornqvist index
      for( i in 1:nGoods) {
         lnp <- lnp + c( 0.5 * ( data[[ shareNames[ i ] ]] +
            baseShares[ i ] *
            matrix( 1, nrow = nObs ) ) * log( data[[ priceNames[ i ] ]] /
            basePrices[ i ] ) )
      }
   } else if(priceIndex=="TL") {      # Translog index
      lnp <- array( coef$alpha0, c( nObs ) )
      for( i in 1:nGoods ) {
         lnp <- lnp + coef$alpha[ i ] * log( data[[ priceNames[ i ] ]] )
         for( j in 1:nGoods ) {
            lnp <- lnp + 0.5 * coef$gamma[ i, j ] *
               log( data[[ priceNames[ i ] ]] ) *
               log( data[[ priceNames[ j ] ]] )
         }
         if( nShifter > 0 ){
            for( j in 1:nShifter ) {
               lnp <- lnp + coef$delta[ i, j ] * data[[ shifterNames[ j ] ]] *
                  log( data[[ priceNames[ i ] ]] )
            }
         }
      }
   } else {
      stop( "the argument 'priceIndex' (price index) must be either 'S',",
         " 'SL', 'P', 'L', 'Ls', 'T' or 'TL'" )
   }

   if( !is.null( row.names( data ) ) ) {
      names( lnp ) <- row.names( data )
   }

   if( exists( "basePrices" ) ){
      attributes( lnp )$basePrices <- basePrices
   }
   if( exists( "baseShares" ) ){
      attributes( lnp )$baseShares <- baseShares
   }

   return( lnp )
}
