cobbDouglasCalc <- function( xNames, data, coef, coefCov = NULL,
      dataLogged = FALSE ) {

   checkNames( c( xNames ), names( data ) )

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

   result <- rep( coef[ "a_0" ], nrow( data ) )
   for( i in seq( along = xNames ) ) {
      result <- result +
         coef[ paste( "a", i, sep = "_" ) ] * logData[[ xNames[ i ] ]]
   }

   if( !dataLogged ) {
      result <- exp( result )
   }

   if( !is.null( coefCov ) ) {
      if( nrow( coefCov ) != nExog + 1 | ncol( coefCov ) != nExog + 1 ) {
         stop( "the covariance matrix of the coefficients",
            " must have exactly ", nExog + 1, " rows and ",
            nExog + 1, " columns" )
      }

      rownames( coefCov ) <- names( coef )
      colnames( coefCov ) <- names( coef )
      jacobian <- matrix( 0, nrow = nrow( data ), ncol = nExog + 1 )
      colnames( jacobian ) <- names( coef )
      jacobian[ , "a_0" ] <- 1
      for( j in 1:nExog ) {
         jacobian[ , paste( "a", j, sep = "_" ) ] <- logData[[ xNames[ j ] ]]
      }
      if( !dataLogged ) {
         jacobian <- jacobian * matrix( rep( result, nExog + 1 ),
            nrow = nrow( data ), ncol = nExog + 1 )
      }
      attributes( result )$variance <-
         diag( jacobian %*% coefCov %*% t( jacobian ) )
   }

   names( result ) <- rownames( data )

   return( result )
}
