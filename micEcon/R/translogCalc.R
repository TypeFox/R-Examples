translogCalc <- function( xNames, data, coef, shifterNames = NULL,
   dataLogged = FALSE ) {

   checkNames( c( xNames, shifterNames ), names( data ) )

   nExog <- length( xNames )
   nShifter <- length( shifterNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2 + nShifter

   if( nCoef > length( coef ) ) {
      stop( "a translog function with ", nExog, " exogenous variables",
         " and ", nShifter, " shifter variables",
         " must have at least ", nCoef, " coefficients" )
   }

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- logDataSet( data = data,
         varNames = xNames, varNamesNum = shifterNames )
   }

   result <- quadFuncCalc( xNames, logData, coef, 
      shifterNames = shifterNames )

   if( !dataLogged ) {
      result <- exp( result )
   }

   return( result )
}
