logDataSet <- function( data, varNames, varNamesNum = NULL ) {

   if( "plm.dim" %in% class( data ) ) {
      logData <- data[ , 1:2 ]
   } else {
      logData <- data.frame( no = c( 1:nrow( data ) ) )
   }
   for( i in seq( along = varNames ) ) {
      logData[[ varNames[ i ] ]] <- log( data[[ varNames[ i ] ]] )
   }
   for( i in seq( along = varNamesNum ) ) {
      if( is.factor( data[[ varNamesNum[ i ] ]] ) | 
            is.logical( data[[ varNamesNum[ i ] ]] ) ) {
         logData[[ varNamesNum[ i ] ]] <- data[[ varNamesNum[ i ] ]]
      } else {
         logData[[ varNamesNum[ i ] ]] <-
            log( data[[ varNamesNum[ i ] ]] )
      }
   }
   if( ! "no" %in% c( varNames, varNamesNum ) ) {
      logData$no <- NULL
   }

   return( logData )
}
