cobbDouglasDeriv <- function( xNames, data, coef, coefCov = NULL,
   yName = NULL, dataLogged = FALSE ) {

   checkNames( c( xNames, yName ), names( data ) )

   nExog <- length( xNames )

   if( nExog + 1 != length( coef ) ) {
      stop( "a Cobb-Douglas function with ", nExog, " exogenous variables",
         " must have exactly ", nExog + 1, " coefficients" )
   }

   coefNames <- paste( "a", c( 0:nExog ), sep = "_" )
   if( is.null( names( coef ) ) ) {
      names( coef ) <- coefNames
   } else {
      coefMissing <- !( coefNames %in% names( coef ) )
      if( any( coefMissing ) ) {
         stop( "coefficient(s) ",
            paste( coefNames[ coefMissing ], collapse = ", " ),
            " are missing" )
      }
      rm( coefMissing )
   }
   rm( coefNames )

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- logDataSet( data = data, varNames = xNames )
   }

   result <- list()

   if( is.null( yName ) ){
      logyHat <- cobbDouglasCalc( xNames = xNames, data = logData,
         coef = coef, dataLogged = TRUE )
   } else {
      if( dataLogged ) {
         logyHat <- data[[ yName ]]
      } else {
         logyHat <- log( data[[ yName ]] )
      }
   }

   deriv <- matrix( NA, nrow( data ), nExog )
   for( i in seq( along = xNames ) ) {
      deriv[ , i ] <- coef[ paste( "a", i, sep = "_" ) ] *
         exp( logyHat ) / exp( logData[[ xNames[ i ] ]] )
   }

   colnames( deriv ) <- xNames
   result$deriv      <- as.data.frame( deriv )

   if( !is.null( coefCov ) ) {
      if( nrow( coefCov ) != nExog + 1 | ncol( coefCov ) != nExog + 1 ) {
         stop( "the covariance matrix of the coefficients",
            " must have exactly ", nExog + 1, " rows and ",
            nExog + 1, " columns" )
      }

      rownames( coefCov ) <- names( coef )
      colnames( coefCov ) <- names( coef )
      variance <- matrix( NA, nrow( data ), nExog )
      for( i in seq( along = xNames ) ) {
         jacobian <- matrix( 0, nrow = nrow( data ), ncol = nExog + 1 )
         colnames( jacobian ) <- names( coef )
         jacobian[ , paste( "a", i, sep = "_" ) ] <-
            exp( logyHat ) / exp( logData[[ xNames[ i ] ]] )
         if( is.null( yName ) ) {
            jacobian[ , "a_0" ] <- coef[ paste( "a", i, sep = "_" ) ] *
               exp( logyHat ) / exp( logData[[ xNames[ i ] ]] )
            for( j in 1:nExog ) {
               jacobian[ , paste( "a", j, sep = "_" ) ] <-
                  jacobian[ , paste( "a", j, sep = "_" ) ] +
                  coef[ paste( "a", i, sep = "_" ) ] * exp( logyHat ) *
                  logData[[ xNames[ j ] ]] / exp( logData[[ xNames[ i ] ]] )
            }
         }
         variance[ , i ] <- diag( jacobian %*% coefCov %*% t( jacobian ) )
      }
      colnames( variance ) <- xNames
      result$variance <- as.data.frame( variance )
   }

   class( result ) <- "cobbDouglasDeriv"
   return( result )
}
