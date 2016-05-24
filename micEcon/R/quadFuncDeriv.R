quadFuncDeriv <- function( xNames, data, coef, coefCov = NULL,
      homWeights = NULL ) {

   # if 'data' is a vector, convert it to a data.frame
   data <- .micEconVectorToDataFrame( data )

   checkNames( c( xNames ), names( data ) )

   # check argument 'homWeights'
   .quadFuncCheckHomWeights( homWeights, xNames )

   result <- list()

   nExog <- length( xNames )

   # check argument 'coef'
   .quadFuncCheckCoefNames( names( coef ), nExog = length( xNames ),
      warn = FALSE )

   # calculate index to normalize variables
   if( !is.null( homWeights ) ) {
      deflator <- 0
      for( i in seq( along = homWeights ) ) {
         deflator <- deflator +
            homWeights[ i ] * data[[ names( homWeights )[ i ] ]]
      }
      whichHom <- which( xNames %in% names( homWeights ) )
   } else {
      whichHom <- NULL
   }

   ## derivatives
   deriv <- array( NA, c( nrow( data ), nExog ) )
   for( i in 1:nExog ) {
      deriv[ , i ] <- coef[ paste( "a", i, sep = "_" ) ]
      for( j in 1:nExog ) {
         deriv[ , i ] <- deriv[ , i ] +
            coef[ paste( "b", min( i, j ), max( i, j ), sep = "_" ) ] * 
            .quadFuncVarHom( data, xNames[ j ], homWeights, deflator )
      }
      if( i %in% whichHom ) {
         deriv[ , i ] <- deriv[ , i ] / deflator
         for( j in whichHom ) {
            deriv[ , i ] <- deriv[ , i ] - homWeights[ xNames[ i ] ] *
               coef[ paste( "a", j, sep = "_" ) ] *
               .quadFuncVarHom( data, xNames[ j ], homWeights, deflator ) / 
               deflator
            for( k in 1:nExog ) {
               deriv[ , i ] <- deriv[ , i ] - homWeights[ xNames[ i ] ] *
                  coef[ paste( "b", min( j, k ), max( j, k ), sep = "_" ) ] *
                  .quadFuncVarHom( data, xNames[ j ], homWeights, deflator ) *
                  .quadFuncVarHom( data, xNames[ k ], homWeights, deflator ) / 
                  deflator
            }
         }
      }
   }
   colnames( deriv ) <- xNames
   deriv    <- as.data.frame( deriv )

   if( !is.null( coefCov ) & is.null( homWeights ) ) {
      ## variances of the derivatives
      variance <- array( NA, c( nrow( data ), nExog ) )
      for(i in 1:nExog ) {
         variance[ , i ] <- coefCov[ paste( "a", i, sep = "_" ), 
            paste( "a", i, sep = "_" ) ]   # variance of aplha(i)
         for( j in 1:nExog ) {
            variance[ , i ] <- variance[ , i ] +
               coefCov[ paste( "a", i, sep = "_" ), 
                  paste( "b", min( i, j ), max( i, j ), sep = "_" ) ] *
               data[[ xNames[ j ] ]]
               # covariance alpha(i)-beta(i,_)
         }
         for( j in 1:nExog ) {
            for( k in 1:nExog ) {
               variance[ , i ] <- variance[ , i ] +
                  coefCov[ paste( "b", min( i, j ), max( i, j ), sep = "_" ),
                     paste( "b", min( i, k ), max( i, k ), sep = "_" ) ] *
                  data[[ xNames[ j ] ]] * data[[ xNames[ k ] ]]
                  # variances + covariance beta(i,_)-beta(i,_)
            }
         }
      }
      stdDev <- variance^0.5  # standard errors
      colnames( variance ) <- xNames
      colnames( stdDev )   <- xNames
      attributes( deriv )$variance <- as.data.frame( variance )
      attributes( deriv )$stdDev   <- as.data.frame( stdDev )
   }

   return( deriv )
}
