translogHessian <- function( xNames, data, coef, yName = NULL,
   dataLogged = FALSE, bordered = FALSE ) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef != length( coef ) ) {
      stop( "a translog function with ", nExog, " exogenous variables",
         " must have exactly ", nCoef, " coefficients" )
   }

   result <- list()

   alpha  <- coef[ 2:( nExog + 1 ) ]
   beta   <- vecli2m( coef[ ( nExog + 2 ):nCoef ] )
   newXNames <- paste( "x.", c( 1:nExog ), sep = "" )
   dNames    <- paste( "d.", c( 1:nExog ), sep = "" )

   if( dataLogged ) {
      logData   <- data.frame( no = c( 1:nrow( data ) ) )
      for( i in seq( along = xNames ) ) {
         logData[[ newXNames[ i ] ]] <- data[[ xNames[ i ] ]]
      }
   } else {
      logData   <- data.frame( no = c( 1:nrow( data ) ) )
      for( i in seq( along = xNames ) ) {
         logData[[ newXNames[ i ] ]] <- log( data[[ xNames[ i ] ]] )
      }
   }

   if( is.null( yName ) ){
      logData$yHat <- translogCalc( newXNames, logData, coef,
         dataLogged = TRUE )
   } else {
      if( dataLogged ) {
         logData$yHat <- data[[ yName ]]
      } else {
         logData$yHat <- log( data[[ yName ]] )
      }
   }

   deriv <- translogDeriv( newXNames, logData, coef, yName = "yHat",
      dataLogged = TRUE )$deriv
   names( deriv ) <- dNames
   logData <- cbind( logData, deriv )

   hessian <- function( values ) {
      result <- matrix( 0, nExog + bordered, nExog + bordered )
      for( i in 1:nExog ) {
         if( bordered ) {
            result[ 1, i + 1 ] <- values[[ dNames[ i ] ]]
            result[ i + 1, 1 ] <- values[[ dNames[ i ] ]]
         }
         for( j in i:nExog ) {
            result[ i + bordered, j + bordered ] <-
               values[[ dNames[ i ] ]] * values[[ dNames[ j ] ]] /
               exp( values[[ "yHat" ]] ) -
               ifelse( i == j, 1, 0 ) * values[[ dNames[ i ] ]] /
               exp( values[[ newXNames[ i ] ]] ) +
               beta[ i, j ] *
               exp( values[[ "yHat" ]] ) /
               ( exp( values[[ newXNames[ i ] ]] ) *
               exp( values[[ newXNames[ j ] ]] ) )
         }
      }
      result[ lower.tri( result ) ] <- t( result )[ lower.tri( result ) ]
      result <- list( result )
      return( result )
   }

   result <- apply( logData, 1, hessian )
   result <- lapply( result, "[[", 1 )
   return( result )
}
