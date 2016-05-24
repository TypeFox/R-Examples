translogRayEst <- function( yNames, xNames, data, shifterNames = NULL, ... ) {

   checkNames( c( yNames, xNames, shifterNames ), names( data ) )

   if( length( yNames ) != 2 ) {
      stop( "the estimation of ray functions has not been implemented",
         " with more than two endogenous variables yet" )
   }

   if( any( c( "distance", "theta" ) %in% c( xNames, shifterNames ) ) ) {
      stop( "the variable names in arguments 'xNames' and 'shifterNames'",
         " must not be 'distance' or 'theta'" )
   }

   nInput <- length( xNames )

   logData <- logDataSet( data = data, varNames = xNames,
      varNamesNum = shifterNames )

   distance <- sqrt( data[[ yNames[ 1 ] ]]^2 + data[[ yNames[ 2 ] ]]^2 )

   logData$distance <- log( distance )

   logData$theta <- acos( data[[ yNames[ 1 ] ]] / distance )

   result <- translogEst( yName = "distance",
      xNames = c( xNames, "theta" ), shifterNames = shifterNames,
      data = logData, dataLogged = TRUE, ... )

   result$call <- match.call()
   result$yName         <- NULL
   result$yNames        <- yNames
   result$xNames        <- xNames
   result$shifterNames  <- shifterNames
   result$distance      <- distance
   result$theta         <- logData$theta

   coefNames <- names( result$coef )
   coefNames[ coefNames ==
      paste( "b", nInput + 1, nInput + 1, sep = "_" ) ] <- "b_t_t"
   coefNames <- sub( paste( "_", nInput + 1, "$", sep = "" ), "_t",
      coefNames )
   names( result$coef ) <- coefNames
   rownames( result$coefCov ) <- coefNames
   colnames( result$coefCov ) <- coefNames

   class( result ) <- c( "translogRayEst", class( result ) )
   return( result )
}
