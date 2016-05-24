translogRayDeriv <- function( yNames, xNames, data, coef,
      dataLogged = FALSE ) {

   checkNames( c( yNames, xNames ), names( data ) )

   if( length( yNames ) > 2 ) {
      stop( "the estimation of ray functions with more than two outputs",
         " has not been implemented yet" )
   } else if( ! length( yNames ) %in% c( 2 ) ) {
      stop( "argument 'yNames' must include exactly two outputs" )
   }

   if( dataLogged ) {
      logData   <- data
   } else {
      logData <- logDataSet( data = data,
         varNames = c( yNames, xNames ) )
   }

   # number of inputs
   nInput <- length( xNames )

   # coefficients with names as if theta was a normal input
   coefTl <- coef
   coefNamesTl <- names( coef )
   coefNamesTl[ coefNamesTl == "b_t_t" ] <-
      paste( "b", nInput + 1, nInput + 1, sep = "_" )
   coefNamesTl <- sub( "_t$", paste( "_", nInput + 1, sep = "" ),
      coefNamesTl )
   names( coefTl ) <- coefNamesTl

   # calculate distance and theta
   distance <- sqrt( exp( logData[[ yNames[ 1 ] ]] )^2 +
         exp( logData[[ yNames[ 2 ] ]] )^2 )
   logData$distance <- log( distance )
   logData$theta <- acos( exp( logData[[ yNames[ 1 ] ]] ) / distance )
   xNamesTheta <- c( xNames, "theta" )

   # elasticities: d log F / d log x_i
   elaData <- translogEla( xNames = xNamesTheta, data = logData, coef = coefTl,
      dataLogged = TRUE )

   # data.frame for derivatives
   deriv <- data.frame( obsNo = rep( NA, nrow( data ) ) )

   # derivatives with respect to the independent variables
   for( i in xNames ) {
      deriv[[ i ]] <- elaData[[ i ]] / exp( logData[[ i ]] )
   }
   if( ! "obsNo" %in% xNames ) {
      deriv$obsNo <- NULL
   }

   # shortcuts for values of y1 and y2
   yVal <- exp( subset( logData, , yNames ) )

   # derivatives with respect to the first dependent variable
   for( j in 1:2 ) {
      deriv[[ yNames[ j ] ]] <- elaData[[ "theta" ]] *
         ( yVal[ , 1 ] * yVal[ , j ] / distance^3 - ( j == 1 ) / distance ) /
         sqrt( 1 - ( yVal[ , 1 ] / distance )^2 ) -
         yVal[ , j ] / distance^2
   }

   rownames( deriv ) <- rownames( data )

   return( deriv )
}