translogDeriv <- function( xNames, data, coef, coefCov = NULL,
   yName = NULL, dataLogged = FALSE ) {

   checkNames( c( xNames, yName ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef != length( coef ) ) {
      stop( "a translog function with ", nExog, " exogenous variables",
         " must have exactly ", nCoef, " coefficients" )
   }

   result <- list()

   alpha  <- coef[ 2:( nExog + 1 ) ]
   beta   <- vecli2m( coef[ ( nExog + 2 ):nCoef ] )

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- logDataSet( data = data,
         varNames = xNames )
   }

   if( is.null( yName ) ){
      logyHat <- translogCalc( xNames, logData, coef,
         dataLogged = TRUE )
   } else {
      if( dataLogged ) {
         logyHat <- data[[ yName ]]
      } else {
         logyHat <- log( data[[ yName ]] )
      }
   }

   deriv <- matrix( 0, nrow( data ), nExog )
   for( i in seq( along = xNames ) ) {
      deriv[ , i ] <- alpha[ i ]
      for( j in seq( along = xNames ) ) {
         deriv[ , i ] <- deriv[ , i ] +
            beta[ i, j ] * logData[[ xNames[ j ] ]]
      }
      deriv[ , i ] <- deriv[ , i ] * exp( logyHat ) / exp( logData[[ xNames[ i ] ]] )
   }

   colnames( deriv ) <- xNames
   result$deriv      <- as.data.frame( deriv )

   class( result ) <- "translogDeriv"
   return( result )
}
