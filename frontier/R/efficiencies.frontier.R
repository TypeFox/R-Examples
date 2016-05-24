# efficiencies of frontier models
efficiencies.frontier <- function( object, asInData = FALSE,
   logDepVar = TRUE, minusU = farrell, farrell = TRUE, 
   margEff = FALSE, ... ) {

   resid <- residuals( object )
   fitted <- - resid
   for( i in 1:nrow( object$dataTable ) ) {
      firm <- object$dataTable[ i, 1 ]
      time <- object$dataTable[ i, 2 ]
      fitted[ firm, time ] <- fitted[ firm, time] + object$dataTable[ i, 3 ]
   }
   sigmaSq <- coef( object )[ "sigmaSq" ]
   gamma <- coef( object )[ "gamma" ]
   lambda <- sqrt( gamma / ( 1 - gamma ) )
   if( minusU ) {
      dir <- 1
   } else {
      dir <- -1
   }
   
   if( object$ineffDecrease ) {
      tau <- 1
   } else {
      tau <- -1
   }

   if( object$modelType == 1 ) {
      if( margEff ) {
         warning( "cannot calculate marginal effects of z variables",
            " for an error components frontier,",
            " because this model does not have z variables" )
         margEff <- FALSE
      }
      if( object$timeEffect ) {
         eta <- coef( object )[ "time" ]
         etaStar <- exp( - eta * ( 1:object$nt - object$nt ) )
      } else {
         eta <- 0
         etaStar <- rep( 1, object$nt )
      }
      if( object$nt == 1 ) {
         residStar <- drop( resid )
         tStar <- 1
      } else {
         residStar <- rep( NA, object$nn )
         tStar <- rep( NA, object$nn )
         for( i in 1:object$nn ) {
            residStar[ i ] <- sum( resid[ i, ] * etaStar, na.rm = TRUE )
            tStar[ i ] <- sum( etaStar[ !is.na( resid[ i, ] ) ]^2 )
         }
      }
      if( object$truncNorm ) {
         mu <- coef( object )[ "mu" ]
      } else {
         mu <- 0
      }
      muStar <- ( - tau * gamma * residStar + mu * ( 1 - gamma ) ) /
         ( 1 + ( tStar - 1 ) * gamma )
      sigmaStarSq <- sigmaSq * gamma * ( 1 - gamma ) /
         ( 1 + ( tStar - 1 ) * gamma )
      sigmaStar <- sqrt( sigmaStarSq )
      result <- matrix( NA, nrow = object$nn,
         ncol = object$nt^object$timeEffect )
      if( logDepVar ) {
         for( j in 1:ncol( result ) ) {
            result[ , j ] <- exp(
               pnorm( - dir * sigmaStar * etaStar[j] + muStar / sigmaStar, 
                  log.p = TRUE ) - pnorm( muStar / sigmaStar, log.p = TRUE ) ) *
               exp( - dir * muStar * etaStar[j] + 0.5 * sigmaStarSq * etaStar[j]^2 )
         }
      } else {
         if( object$ineffDecrease == minusU ) {
            fittedStar <- rep( NA, object$nn )
            tInd <- rep( NA, object$nn )
            for( i in 1:object$nn ) {
               fittedStar[ i ] <- sum( fitted[ i, ], na.rm = TRUE )
               tInd[ i ] <- sum( !is.na( resid[ i, ] ) )
            }
            for( j in 1:ncol( result ) ) {
               result[ , j ] <- 1 - tau * etaStar[ j ] * ( muStar + sigmaStar *
                  exp( dnorm( muStar / sigmaStar, log = TRUE ) -
                     pnorm( muStar / sigmaStar, log.p = TRUE ) ) ) /
                  ( fittedStar / tInd )
            }
         } else {
            warning( "currently, the efficiency estimates based on models",
               " with non-logged dependent variable can be calculated only",
               " if 'ineffDecrease' is equal to 'minusU'" )
         }
      }
      # set efficiency estimates of missing observations to NA
      if( ncol( result ) > 1 ) {
         result[ is.na( resid ) ] <- NA
      }
   } else if( object$modelType == 2 ) {
      if( object$zIntercept ) {
         zDelta <- coef( object )[ "Z_(Intercept)" ]
      } else {
         zDelta <- 0
      }
      nz <- ncol( object$dataTable ) - 3 - object$nb
      if( nz > 0 ) {
         for( i in 1:nz ) {
            zDelta <- zDelta + object$dataTable[ ,
                  ncol( object$dataTable ) - nz + i ] *
               coef( object )[ object$icept + object$nb + object$zIntercept + i  ]
         }
      } else {
         zDelta <- rep( zDelta, nrow( object$dataTable ) )
      }
      zDeltaMat <- matrix( NA, nrow = object$nn, ncol = object$nt )
      for( i in 1:object$nob ) {
         zDeltaMat[ object$dataTable[ i, 1 ], object$dataTable[ i, 2 ] ] <-
            zDelta[ i ]
      }
      sigmaBarSq <- gamma * ( 1 - gamma ) * sigmaSq
      sigmaBar <- sqrt( sigmaBarSq )
      muBar <- ( 1 - gamma ) * zDeltaMat - tau * gamma * resid
      if( logDepVar ) {
         result <- exp( pnorm( - dir * sigmaBar + muBar / sigmaBar, log.p = TRUE ) -
               pnorm( muBar / sigmaBar, log.p = TRUE ) ) *
               exp( - dir * muBar + 0.5 * sigmaBarSq )
            
         if( margEff ) {
            if( nz < 1 ) {
               warning( "cannot calculate marginal effects of z variables",
                  " for a model that does not have z variables" )
               margEff <- FALSE
            } else {
               margEffectsBase <- ( ( exp( 
                     dnorm( - dir * sigmaBar + muBar / sigmaBar, log = TRUE ) -
                     pnorm( muBar / sigmaBar, log.p = TRUE ) ) / sigmaBar ) *
                  exp( - dir * muBar + 0.5 * sigmaBarSq ) +
                  exp( pnorm( - dir * sigmaBar + muBar / sigmaBar, log.p = TRUE ) -
                     pnorm( muBar / sigmaBar, log.p = TRUE ) ) *
                  exp( - dir * muBar + 0.5 * sigmaBarSq ) * ( - dir ) -
                  exp( dnorm( muBar / sigmaBar, log = TRUE ) +
                     pnorm( - dir * sigmaBar + muBar / sigmaBar, log.p = TRUE ) -
                     2 * pnorm( muBar / sigmaBar, log.p = TRUE ) ) *
                  ( exp( - dir * muBar + 0.5 * sigmaBarSq ) / sigmaBar ) ) *
                  ( 1 - gamma )
               margEffects <- array( NA, 
                  c( nrow( margEffectsBase ), ncol( margEffectsBase ), nz ) )
               for( i in 1:nz ) {
                  margEffects[ , , i ] <- margEffectsBase *
                     coef( object )[ object$icept + object$nb + object$zIntercept + i  ]
               }
            }
         }
      } else {
         if( object$ineffDecrease == minusU ) {
            result <- 1 - tau * ( muBar + sigmaBar *
               exp( dnorm( muBar / sigmaBar, log = TRUE ) -
                  pnorm( muBar / sigmaBar, log.p = TRUE ) ) ) /
               fitted
         } else {
            result <- matrix( NA, nrow = nrow( fitted ), ncol = ncol( fitted ) )
            warning( "currently, the efficiency estimates based on models",
               " with non-logged dependent variable can be calculated only",
               " if 'ineffDecrease' is equal to 'minusU'" )
         }
         if( margEff ) {
            warning( "calculation of marginal effects of z variables",
               " has not been implemented for models with non-logged",
               " dependent variables yet" )
            margEff <- FALSE
         }
      }
   } else {
      stop( "internal error: unknow model type '",
         object$modelType, "'" )
   }

   if( minusU ) {
      result[ result > 1 ] <- 1
   } else {
      result[ result < 1 ] <- 1
   }

   rownames( result ) <- rownames( resid )
   if( ncol( result ) > 1 ) {
      colnames( result ) <- colnames( resid )
   } else {
      colnames( result ) <- "efficiency"
   }

   if( margEff ) {
      dimnames( margEffects ) <- list( rownames( resid ),
         if( ncol( result ) > 1 ){ colnames( resid ) } else { "efficiency" },
         names( coef( object ) )[ ( object$icept + object$nb + object$zIntercept + 1 ):( 
            object$icept + object$nb + object$zIntercept + nz ) ] )
   }

   if( asInData ) {
      effic <- rep( NA, length( object$validObs ) )
      if( sum( object$validObs ) != nrow( object$dataTable ) ) {
         stop( "internal error: number of rows of element 'dataTable' is not",
            " equal to number of valid observations" )
      }
      for( i in 1:nrow( object$dataTable ) ) {
         effic[ object$validObs ][ i ] <- result[ object$dataTable[ i , 1 ],
            min( object$dataTable[ i , 2 ], ncol( result ) ) ]
      }
      result <- effic
      names( result ) <- names( object$validObs )
      if( margEff ) {
         margEffects2 <- matrix( NA, nrow = length( object$validObs ), ncol = nz )
         for( i in 1:nrow( object$dataTable ) ) {
            for( k in 1:nz ) { 
               margEffects2[ object$validObs, k ][ i ] <- 
                  margEffects[ object$dataTable[ i , 1 ],
                  min( object$dataTable[ i , 2 ], ncol( result ) ), k ]
            }
         }
         rownames( margEffects2 ) <- names( object$validObs )
         colnames( margEffects2 ) <- dimnames( margEffects )[[ 3 ]]
         margEffects <- margEffects2
      }
   }

   if( margEff ) {
      attr( result, "margEff" ) <- margEffects
   }

   return( result )
}
