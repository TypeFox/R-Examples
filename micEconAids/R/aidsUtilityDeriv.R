aidsUtilityDeriv <- function( priceNames, totExpName, coef, data,
      rel = FALSE ) {

   # check argument 'coef' (coefficients)
   coefCheckResult <- .aidsCheckCoef( coef, variables = list(
      list( length( priceNames ), "prices", "goods"  ) ) )
   if( !is.null( coefCheckResult ) ){
      stop( coefCheckResult )
   }
   if( is.null( coef$alpha0 ) ) {
      stop( "argument 'coef' has no element 'alpha0'" )
   }
   if( is.null( coef$beta0 ) ) {
      coef$beta0 <- 1
   }

   # checking argument 'data'
   if( class( data ) != "data.frame" ) {
      stop( "argument 'data' must be a data frame" )
   }

   # check names of variables
   checkNames( c( priceNames, totExpName ), names( data ) )

   # number of goods
   nGoods <- length( priceNames )

   # data.frame that contains the partial derivatives
   result <- data.frame( temp = rep( NA, nrow( data ) ) )

   termA <- rep( coef$alpha0, nrow( result ) )
   for( i in 1:nGoods ) {
      termA <- termA + coef$alpha[ i ] * log( data[[ priceNames[ i ] ]] )
      for( j in 1:nGoods ) {
         termA <- termA + 0.5 * coef$gamma[ i, j ] *
            log( data[[ priceNames[ i ] ]] ) *
            log( data[[ priceNames[ j ] ]] )
      }
   }

   termB <- rep( coef$beta0, nrow( result ) )
   for( i in 1:nGoods ) {
      termB <- termB * data[[ priceNames[ i ] ]]^coef$beta[ i ]
   }

   if( rel ) {
      utility <- aidsUtility( priceNames = priceNames,
         totExpName = totExpName, coef = coef, data = data )
   }

   for( i in 1:nGoods ){
      numerator <- - coef$alpha[ i ] -
         coef$beta[ i ] * ( log( data[[ totExpName ]] ) - termA )
      for( j in 1:nGoods ) {
         numerator <- numerator -
            0.5 * ( coef$gamma[ i, j ] + coef$gamma[ j, i ] ) *
            log( data[[ priceNames[ j ] ]] )
      }
      denominator <- data[[ priceNames[ i ] ]] * termB
      result[[ priceNames[ i ] ]] <- numerator / denominator
      if( rel ) {
         result[[ priceNames[ i ] ]] <- result[[ priceNames[ i ] ]] *
            data[[ priceNames[ i ] ]] / utility
      }
   }

   result[[ totExpName ]] <- 1 / ( data[[ totExpName ]] * termB )
   if( rel ) {
      result[[ totExpName ]] <- result[[ totExpName ]] *
         data[[ totExpName ]] / utility
   }

   result$temp <- NULL

   return( result )
}