frontierTranslogRay <- function( yNames, xNames, shifterNames = NULL,
      zNames = NULL, data, ... ) {

   checkNames( c( yNames, xNames, shifterNames, zNames ), names( data ) )

   nOutput <- length( yNames )

   if( nOutput < 2 ) {
      stop( "the argument 'yNames' must include the names of",
         "at least two output variables" )
   }

   if( any( c( "distance", "theta" ) %in% c( xNames, shifterNames ) ) ) {
      stop( "the variable names in arguments 'xNames' and 'shifterNames'",
         " must not be 'distance' or 'theta'" )
   }

   nInput <- length( xNames )

   logData <- logDataSet( data = data, varNames = xNames )

   distance <- 0
   for( i in 1:nOutput ) {
      distance <- distance + data[[ yNames[ i ] ]]^2
   }
   distance <- sqrt( distance )

   logData$distance <- log( distance )

   sinProd <- 1
   for( i in 1:( nOutput - 1 ) ) {
      ratio <- data[[ yNames[ i ] ]] / ( distance * sinProd )
      ratio[ ratio > 1 ] <- 1
      ratio[ ratio < -1 ] <- -1
      logData[[ paste( "theta", i, sep = "_" ) ]] <- acos( ratio )
      sinProd <- sinProd * sin( logData[[ paste( "theta", i, sep = "_" ) ]] )
   }

   # shifter variables
   for( i in seq( along = shifterNames ) ) {
      logData[[ shifterNames[ i ] ]] <- data[[ shifterNames[ i ] ]]
   }

   # z variables
   for( i in seq( along = zNames ) ) {
      logData[[ zNames[ i ] ]] <- data[[ zNames[ i ] ]]
   }

   result <- frontierQuad( yName = "distance",
      xNames = c( xNames, paste( "theta", 1:( nOutput - 1 ), sep = "_" ) ), 
      shifterNames = shifterNames,
      zNames = zNames, data = logData, ... )

   result$call <- match.call()
   result$yName         <- NULL
   result$yNames        <- yNames
   result$xNames        <- xNames
   result$shifterNames  <- shifterNames
   result$distance      <- distance
   for( i in 1:( nOutput - 1 ) ) {
      result[[ paste( "theta", i, sep = "_" ) ]] <- 
         logData[[ paste( "theta", i, sep = "_" ) ]]
   }

   coefNames <- names( result$mleParam )[
      1:( 1 + ( nInput + nOutput - 1 ) + 
         ( nInput + nOutput ) * ( nInput + nOutput - 1 ) / 2 ) ]

   for( i in 1:( nOutput - 1 ) ) {
      coefNames <- gsub( paste( "_", nInput + i, "$", sep = "" ), 
         paste( "_t", i, sep = "" ), coefNames )
      coefNames <- gsub( paste( "_", nInput + i, "_", sep = "" ), 
         paste( "_t", i, "_", sep = "" ), coefNames )
   }
   names( result$olsParam )[ 1:length( coefNames ) ] <- coefNames
   names( result$olsStdEr )[ 1:length( coefNames ) ] <- coefNames
   if( ! is.null( result$gridParam ) ) {
      names( result$gridParam )[ 1:length( coefNames ) ] <- coefNames
   }
   names( result$mleParam )[ 1:length( coefNames ) ] <- coefNames
   rownames( result$mleCov )[ 1:length( coefNames ) ] <- coefNames
   colnames( result$mleCov )[ 1:length( coefNames ) ] <- coefNames

   class( result ) <- c( "frontierTranslogRay", class( result ) )

   return( result )
}
