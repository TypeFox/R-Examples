quadFuncEst <- function( yName, xNames, data, shifterNames = NULL,
   linear = FALSE, homWeights = NULL, regScale = 1, ... ) {

   checkNames( c( yName, xNames, shifterNames ), names( data ) )

   # check argument 'homWeights'
   .quadFuncCheckHomWeights( homWeights, xNames )

   nExog   <- length( xNames )
   nShifter <- length( shifterNames )
   result <- list()
   result$call <- match.call()

   isPanel <- any( c( "plm.dim", "pdata.frame" ) %in% class( data ) )

   if( isPanel ) {
      estData <- data[ , 1:2 ]
      estData$y <- data[[ yName ]]
   } else {
      estData <- data.frame( y = data[[  yName ]] )
   }

   if( !is.null( homWeights ) ) {
      estData$deflator <- 0
      for( i in seq( along = homWeights ) ) {
         estData$deflator <- estData$deflator + 
            homWeights[ i ] * data[[ names( homWeights )[ i ] ]]
      }
      whichHom <- which( xNames %in% names( homWeights ) )
      xOmit <- names( homWeights )[ 1 ]
      iOmit <- which( xNames == xOmit )
   } else {
      iOmit <- 0
      xOmit <- NULL
   }

   estFormula <- "y ~ 1"
   for( i in seq( along = xNames ) ) {
      if( i != iOmit ) {
         xName <- paste( "a", as.character( i ), sep = "_" )
         estData[[ xName ]] <- .quadFuncVarHom( data, xNames[ i ], 
            homWeights, estData$deflator, xOmit ) / regScale
         estFormula <- paste( estFormula, "+", xName )
      }
   }
   if( !linear ) {
      for( i in seq( along = xNames ) ) {
         for( j in i:nExog ) {
            if( i != iOmit & j != iOmit ) {
               xName <- paste( "b", as.character( i ), as.character( j ),
                  sep = "_" )
               estData[[ xName ]] <- 0.5 *
                  ifelse( i == j, 1, 2 ) *
                  .quadFuncVarHom( data, xNames[ i ], homWeights, 
                     estData$deflator, xOmit ) * 
                  .quadFuncVarHom( data, xNames[ j ], homWeights, 
                     estData$deflator, xOmit ) / 
                  regScale
               estFormula <- paste( estFormula, "+", xName )
            }
         }
      }
   }
   for( i in seq( along = shifterNames ) ) {
      if( is.factor( data[[ shifterNames[ i ] ]] ) | 
            is.logical( data[[ shifterNames[ i ] ]] ) ) {
         xName <- paste( "d", "_", as.character( i ), "_", sep = "" )
         estData[[ xName ]] <- data[[ shifterNames[ i ] ]]
      } else {
         xName <- paste( "d", as.character( i ), sep = "_" )
         estData[[ xName ]] <- data[[ shifterNames[ i ] ]] / regScale
      }
      estFormula <- paste( estFormula, "+", xName )
   }
   result$nExog <- nExog
   result$nShifter <- nShifter
   if( isPanel ) {
      result$est <- plm( as.formula( estFormula ), estData, ... )
      result$est$call$formula <- as.formula( estFormula )
   } else {
      result$est <- lm( as.formula( estFormula ), estData, ... )
   }
   result$residuals <- c( residuals( result$est ) )
   result$fitted    <- estData$y - result$residuals

   # coefficients and their covariance matrix
   result$coef      <- coef( result$est )
   result$coefCov   <- vcov( result$est )
   if( isPanel ) {
      if( is.null( result$est$call$model ) ||
            result$est$call$model == "within" ) {
         result$coef <- c( mean( fixef( result$est ) ), result$coef )
         result$coefCov <- rbind( NA, cbind( NA, vcov( result$est ) ) )
      }
   }
   names( result$coef )[ 1 ]       <- "a_0"
   rownames( result$coefCov )[ 1 ] <- "a_0"
   colnames( result$coefCov )[ 1 ] <- "a_0"

   # adding coefficients and covariances that have been dropped 
   # due to the homogeneity restriction
   if( !is.null( homWeights ) ) {
      # missing coefficients
      coefOmit <- 0
      for( i in whichHom[ whichHom != iOmit ] ) {
         coefOmit <- coefOmit - result$coef[ paste( "a", i, sep = "_" ) ]
      }
      result$coef <- c( result$coef, coefOmit )
      names( result$coef )[ length( result$coef ) ] <- 
         paste( "a", iOmit, sep = "_" )
      if( !linear & nExog > 0 ) {
         for( i in c( (1:nExog)[ (1:nExog) != iOmit ], iOmit ) ) {
            coefOmit <- 0
            for( j in whichHom[ whichHom != iOmit ] ) {
               coefOmit <- coefOmit - result$coef[ 
                  paste( "b", min( i, j ), max( i, j ), sep = "_" ) ]
            }
            result$coef <- c( result$coef, coefOmit )
            names( result$coef )[ length( result$coef ) ] <- 
               paste( "b", min( i, iOmit ), max( i, iOmit ), sep = "_" )
         }
      }
      # missing rows of covariance matrix
      coefCovOmit <- rep( 0, ncol( result$coefCov ) )
      for( i in whichHom[ whichHom != iOmit ] ) {
         coefCovOmit <- coefCovOmit - 
            result$coefCov[ paste( "a", i, sep = "_" ), ]
      }
      result$coefCov <- rbind( result$coefCov, coefCovOmit )
      rownames( result$coefCov )[ nrow( result$coefCov ) ] <- 
         paste( "a", iOmit, sep = "_" )
      if( !linear & nExog > 0 ) {
         for( i in c( (1:nExog)[ (1:nExog) != iOmit ], iOmit ) ) {
            coefCovOmit <- rep( 0, ncol( result$coefCov ) )
            for( j in whichHom[ whichHom != iOmit ] ) {
               coefCovOmit <- coefCovOmit - result$coefCov[ 
                  paste( "b", min( i, j ), max( i, j ), sep = "_" ), ]
            }
            result$coefCov <- rbind( result$coefCov, coefCovOmit )
            rownames( result$coefCov )[ nrow( result$coefCov ) ] <- 
               paste( "b", min( i, iOmit ), max( i, iOmit ), sep = "_" )
         }
      }
      # missing columns of covariance matrix
      coefCovOmit <- rep( 0, nrow( result$coefCov ) )
      for( i in whichHom[ whichHom != iOmit ] ) {
         coefCovOmit <- coefCovOmit - 
            result$coefCov[ , paste( "a", i, sep = "_" ) ]
      }
      result$coefCov <- cbind( result$coefCov, coefCovOmit )
      colnames( result$coefCov )[ ncol( result$coefCov ) ] <- 
         paste( "a", iOmit, sep = "_" )
      if( !linear & nExog > 0 ) {
         for( i in c( (1:nExog)[ (1:nExog) != iOmit ], iOmit ) ) {
            coefCovOmit <- rep( 0, nrow( result$coefCov ) )
            for( j in whichHom[ whichHom != iOmit ] ) {
               coefCovOmit <- coefCovOmit - result$coefCov[ ,
                  paste( "b", min( i, j ), max( i, j ), sep = "_" ) ]
            }
            result$coefCov <- cbind( result$coefCov, coefCovOmit )
            colnames( result$coefCov )[ ncol( result$coefCov ) ] <- 
               paste( "b", min( i, iOmit ), max( i, iOmit ), sep = "_" )
         }
      }
   }

   if( linear & nExog > 0 ) {
      nQuadCoef <- nExog * ( nExog + 1 ) / 2
      quadCoefNames <- paste( "b", 
         vecli( matrix( rep( 1:nExog, nExog ), nrow = nExog ) ), 
         vecli( matrix( rep( 1:nExog, each = nExog ), nrow = nExog ) ),
         sep = "_" )
      quadCoef <- rep( 0, nQuadCoef )
      names( quadCoef ) <- quadCoefNames
      result$coef <- c( result$coef, quadCoef )
      quadCoefCovRows <- matrix( 0, nrow = nQuadCoef, 
         ncol = ncol( result$coefCov ) ) 
      rownames( quadCoefCovRows ) <- quadCoefNames
      result$coefCov <- rbind( result$coefCov, quadCoefCovRows )
      quadCoefCovCols <- matrix( 0, nrow = nrow( result$coefCov ), 
         ncol = nQuadCoef ) 
      colnames( quadCoefCovCols ) <- quadCoefNames
      result$coefCov <- cbind( result$coefCov, quadCoefCovCols )
   }

   result$coef <- result$coef[ .micEconCoefOrder( names( result$coef ) ) ]
   result$coefCov <- result$coefCov[
      .micEconCoefOrder( rownames( result$coefCov ) ), ]
   result$coefCov <- result$coefCov[ ,
      .micEconCoefOrder( colnames( result$coefCov ) ) ]

   if( isPanel ) {
      result$r2    <- unname( summary( result$est )$r.squared[ "rsq" ] )
      result$r2bar <- unname( summary( result$est )$r.squared[ "adjrsq" ] )
   } else {
      result$r2    <- summary( result$est )$r.squared
      result$r2bar <- summary( result$est )$adj.r.squared
   }
   result$nObs  <- length( result$residuals )
   result$yName        <- yName
   result$xNames       <- xNames
   result$shifterNames <- shifterNames
   result$homWeights   <- homWeights
   result$regScale     <- regScale

   if( !isPanel ) {
      result$model.matrix <- model.matrix( result$est )
   }
   class( result ) <- "quadFuncEst"
   return( result )
}
