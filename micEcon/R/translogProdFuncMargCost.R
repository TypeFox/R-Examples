translogProdFuncMargCost <- function( yNames, xNames, wNames,
      data, coef, dataLogged = FALSE ) {

   checkNames( c( yNames, xNames, wNames ), names( data ) )

   if( length( yNames ) > 2 ) {
      stop( "the estimation of ray functions has not been implemented",
         " with more than two outputs yet" )
   } else if( ! length( yNames ) %in% c( 1, 2 ) ) {
      stop( "argument 'yNames' must include either one or two outputs" )
   }

   if( length( xNames ) != length( wNames ) ) {
      stop( "arguments 'xNames' and 'wNames' must have the same length" )
   }

   if( dataLogged ) {
      logData   <- data
   } else {
      logData <- logDataSet( data = data,
         varNames = c( yNames, xNames, wNames ) )
   }

   # number of inputs
   nInput <- length( xNames )

   # coefficients with names as if theta was a normal input
   coefTl <- coef
   if( length( yNames ) > 1 ) {
      coefNamesTl <- names( coef )
      coefNamesTl[ coefNamesTl == "b_t_t" ] <-
         paste( "b", nInput + 1, nInput + 1, sep = "_" )
      coefNamesTl <- sub( "_t$", paste( "_", nInput + 1, sep = "" ),
         coefNamesTl )
      names( coefTl ) <- coefNamesTl
   }

   # calculate distance and theta
   xNamesTheta <- xNames
   if( length( yNames ) >  1 ) {
      distance <- sqrt( exp( logData[[ yNames[ 1 ] ]] )^2 +
            exp( logData[[ yNames[ 2 ] ]] )^2 )
      logData$distance <- log( distance )
      logData$theta <- acos( exp( logData[[ yNames[ 1 ] ]] ) / distance )
      xNamesTheta <- c( xNamesTheta, "theta" )
   }

   # elasticities: d log F / d log x_i
   elaData <- translogEla( xNames = xNamesTheta, data = logData, coef = coefTl,
      dataLogged = TRUE )

   # matrix of beta coefficients
   beta <- matrix( NA, nrow = nInput, ncol = nInput )
   for( i in 1:nInput ) {
      for( j in i:nInput ) {
         beta[ i, j ] <- coef[ paste( "b", i, j, sep = "_" ) ]
         beta[ j, i ] <- beta[ i, j ]
      }
   }

   # marginal costs
   margCost <- matrix( NA, nrow = nrow( logData ), ncol = length( yNames ) )

   for( t in 1:nrow( logData ) ) {
      yVal <- exp( as.numeric( logData[ t, yNames ] ) )
      xVal <- exp( as.numeric( logData[ t, xNames ] ) )
      wVal <- exp( as.numeric( logData[ t, wNames ] ) )
      eVal <- as.numeric( elaData[ t, ] )
      if( length( yNames ) > 1 ) {
         dVal <- distance[ t ]
      }

      # Jacobian matrix of g with respect to x
      gxJac <- matrix( NA, nrow = nInput, ncol = nInput )
      for( j in 1:nInput ) {
         for( i in 1:( nInput - 1 ) ) {
            gxJac[ i, j ] <-
               ( i == j ) * wVal[ nInput ] *
                  ( xVal[ nInput ] / xVal[ j ]^2 ) *
                  eVal[ i ] / eVal[ nInput ] -
               ( j == nInput ) * wVal[ nInput ] * ( 1 / xVal[ i ] ) *
                  eVal[ i ] / eVal[ nInput ] -
               wVal[ nInput ] * ( xVal[ nInput ] / xVal[ i ] ) *
                  ( beta[ i, j ] / xVal[ j ] ) / eVal[ nInput ] +
               wVal[ nInput ] * ( xVal[ nInput ] / xVal[ i ] ) *
                  ( eVal[ i ] / eVal[ nInput ]^2 ) *
                  beta[ nInput, j ] / xVal[ j ]
         }
         gxJac[ nInput, j ] <- eVal[ j ] / xVal[ j ]
      }

      # Jacobian matrix of g with respect to y
      gyJac <- matrix( nrow = nInput, ncol = length( yNames ) )
      if( length( yNames ) == 1 ) {
         gyJac[ 1:( nInput - 1 ), 1 ] <- 0
         gyJac[ nInput, 1 ] <- - 1 / yVal
      } else {
         # derivatives of theta with respect to the outputs
         derivThetaY <- rep( NA, length( yNames ) )
         for( j in 1:length( yNames ) ) {
            derivThetaY[ j ] <-
               ( yVal[ 1 ] * yVal[ j ] / dVal^3 - ( j == 1 ) / dVal ) /
               ( sqrt( 1 - ( yVal[ 1 ] / dVal )^2 ) )
         }
         # derivatives of theta with respect to the outputs
         derivDistanceY <- rep( NA, length( yNames ) )
         for( j in 1:length( yNames ) ) {
            derivDistanceY[ j ] <- yVal[ j ] / dVal^2
         }
         # and now the Jacobian ...
         for( j in 1:length( yNames ) ) {
            for( i in 1:( nInput - 1 ) ) {
               gyJac[ i, j ] <-
                  ( derivThetaY[ j ] * wVal[ nInput ] * xVal[ nInput ] /
                     ( xVal[ i ] * eVal[ nInput ] ) ) *
                  ( ( eVal[ i ] / eVal[ nInput ] ) *
                     coef[ paste( "b", nInput, "t", sep = "_" ) ] - 
                     coef[ paste( "b", i, "t", sep = "_" ) ] )
            }
            gyJac[ nInput, j ] <- eVal[ nInput + 1 ] * derivThetaY[ j ] -
               - derivDistanceY[ j ]
         }
      }

      # Jacobian matrix of x with respect to y
      xyJac <- - solve( gxJac, gyJac )

      # marginal costs
      margCost[ t, ] <- t( wVal ) %*% xyJac
   }

   margCost <- as.data.frame( margCost )
   rownames( margCost ) <- rownames( data )
   names( margCost ) <- yNames

   return( margCost )
}