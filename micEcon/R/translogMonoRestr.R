translogMonoRestr <- function( xNames, data,
   dataLogged = FALSE, box = FALSE ) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- logDataSet( data = data,
         varNames = xNames )
   }

   if( box ) {
      extremeLogValues <- list()
      for( i in seq( along = xNames ) ) {
         extremeLogValues[[ i ]] <- c(
            min( logData[[ xNames[ i ] ]] ),
            max( logData[[ xNames[ i ] ]] ) )
      }
      logData <- expand.grid( extremeLogValues )
      colnames( logData ) <- xNames
   }
   nObs  <- nrow( logData )
   restr <- matrix( 0, nObs * nExog, nCoef )
   for( i in seq( along = xNames ) ) {
      myRows <- c( ( ( i - 1 ) * nObs + 1 ):( i * nObs ) )
      restr[ myRows, 1 + i ] <- 1
      for( j in seq( along = xNames ) ) {
         restr[ myRows, 1 + nExog + veclipos( i, j, nExog ) ] <-
            logData[[ xNames[ j ] ]]
      }
   }

   return( restr )
}

# test with (only if box == FALSE):
# matrix( translogMonoRestr( lnInputNames, estData ) %*% coef( a ), ncol=4) *
# exp(translogCalc( lnInputNames, estData, coef( a ) ) %*% t(c(1,1,1,1))) /
# estData[ , inputNames ] -
# translogDeriv( lnInputNames, estData, coef( a ) )$deriv
