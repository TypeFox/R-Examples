aidsUtility <- function( priceNames, totExpName, coef, data ) {

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

   numerator <- log( data[[ totExpName ]] ) - coef$alpha0
   for( i in 1:nGoods ) {
      numerator <- numerator -
         coef$alpha[ i ] * log( data[[ priceNames[ i ] ]] )
      for( j in 1:nGoods ) {
         numerator <- numerator - 0.5 * coef$gamma[ i, j ] *
            log( data[[ priceNames[ i ] ]] ) *
            log( data[[ priceNames[ j ] ]] )
      }
   }

   denominator <- rep( coef$beta0, length( numerator ) )
   for( i in 1:nGoods ) {
      denominator <- denominator * data[[ priceNames[ i ] ]]^coef$beta[ i ]
   }

   result <- numerator / denominator

   return( result )
}