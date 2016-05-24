invMillsRatio <- function( x, all = FALSE ) {
   errorMessage <- paste( "calculating the 'Inverse Mills Ratio' only works",
      "for probit models estimated by 'glm' or 'vglm' or 'probit'." )
   if(inherits( x, "glm") ) {
      if( x$family$family != "binomial" || x$family$link != "probit" ) {
         stop( errorMessage )
      }
      result <- data.frame( no = 1:nrow( model.frame( x ) ),
         row.names = rownames( model.frame( x ) ) )
      result$IMR1 <- dnorm( x$linear.predictors ) /
         pnorm( x$linear.predictors )
      result$delta1 <- result$IMR1 * ( result$IMR1 + x$linear.predictors )
      result$IMR0 <- dnorm( x$linear.predictors ) /
         pnorm( -x$linear.predictors )
      result$delta0 <- result$IMR0 * ( result$IMR0 + x$linear.predictors )
   } else if(inherits(x, "probit")) {
      # Note: 'probit' need not to be the first component in the class
      result <- data.frame( no = seq(length=nObs(x)),
         row.names = rownames( model.matrix( x ) ) )
      result$IMR1 <- dnorm(linearPredictors(x))/pnorm(linearPredictors(x))
      result$delta1 <- result$IMR1 * ( result$IMR1 + linearPredictors(x))
      result$IMR0 <- dnorm(linearPredictors(x))/pnorm(-linearPredictors(x))
      result$delta0 <- result$IMR0 * ( result$IMR0 + linearPredictors(x))

   } else if(inherits(x, "vglm")) {
      if( x@family@blurb[1] != "Bivariate probit model\n"  ) {
         stop( errorMessage )
      }
      result <- data.frame( no = 1:nrow( x@predictors ),
         row.names = rownames( x@predictors ) )
      if( length( x@misc$link ) == 1 ) {
         vglmLink <- x@misc$link
      } else {
         vglmLink <- x@misc$link[ "rho" ]
      }
      if( vglmLink %in% c( "identity", "identitylink" ) ) {
         rho <- x@predictors[ , 3 ]
      } else if( vglmLink == "rhobit" ){
         rho <- rhobit( x@predictors[ , 3 ], inverse = TRUE )
      } else {
         stop( "the bivariate probit (binom2.rho) must be either estimated",
            " with link 'rhobit' or 'identity'/'identitylink'" )
      }
      if( max( rho ) > 1 ) {
         stop( "the correlation between the error terms (rho) is larger",
            " than 1" )
      }
      pmvnormValues11 <- rep( NA, nrow( result ) )
      pmvnormValues10 <- rep( NA, nrow( result ) )
      pmvnormValues01 <- rep( NA, nrow( result ) )
      pmvnormValues00 <- rep( NA, nrow( result ) )
      for( i in 1:nrow( result ) ) {
         corrEq <- matrix( c( 1, rho[ i ], rho[ i ], 1 ), ncol = 2 )
         corrUneq <- matrix( c( 1, -rho[ i ], -rho[ i ], 1 ), ncol = 2 )
         if( x@y[ i, "11" ] == 1 | all ) {
            pmvnormValues11[ i ] <- pmvnorm( upper = x@predictors[ i, 1:2 ],
               corr = corrEq )
         }
         if( x@y[ i, "10" ] == 1 | all ) {
            pmvnormValues10[ i ] <- pmvnorm( upper = c( x@predictors[ i, 1 ],
               -x@predictors[ i, 2 ] ), corr = corrUneq )
         }
         if( x@y[ i, "01" ] == 1 | all ) {
            pmvnormValues01[ i ] <- pmvnorm( upper = c( -x@predictors[ i, 1 ],
               x@predictors[ i, 2 ] ), corr = corrUneq )
         }
         if( x@y[ i, "00" ] == 1 | all ) {
            pmvnormValues00[ i ] <- pmvnorm( upper = -x@predictors[ i, 1:2 ],
               corr = corrEq )
         }
      }

      result$IMR11a <- dnorm( x@predictors[ , 1 ] ) *
         pnorm( ( x@predictors[ , 2 ] - rho * x@predictors[ , 1 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues11
      result$IMR11b <- dnorm( x@predictors[ , 2 ] ) *
         pnorm( ( x@predictors[ , 1 ] - rho * x@predictors[ , 2 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues11
      result$IMR10a <- dnorm( x@predictors[ , 1 ] ) *
         pnorm( ( -x@predictors[ , 2 ] + rho * x@predictors[ , 1 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues10
      result$IMR10b <- -dnorm( x@predictors[ , 2 ] ) *
         pnorm( ( x@predictors[ , 1 ] - rho * x@predictors[ , 2 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues10
      result$IMR01a <- -dnorm( x@predictors[ , 1 ] ) *
         pnorm( ( x@predictors[ , 2 ] - rho * x@predictors[ , 1 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues01
      result$IMR01b <- dnorm( x@predictors[ , 2 ] ) *
         pnorm( ( -x@predictors[ , 1 ] + rho * x@predictors[ , 2 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues01
      result$IMR00a <- -dnorm( x@predictors[ , 1 ] ) *
         pnorm( ( -x@predictors[ , 2 ] + rho * x@predictors[ , 1 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues00
      result$IMR00b <- -dnorm( x@predictors[ , 2 ] ) *
         pnorm( ( -x@predictors[ , 1 ] + rho * x@predictors[ , 2 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues00

      if( all ) {
         selection1X <- rep( TRUE, nrow( result ) )
         selection0X <- rep( TRUE, nrow( result ) )
         selectionX1 <- rep( TRUE, nrow( result ) )
         selectionX0 <- rep( TRUE, nrow( result ) )
      } else {
         selection1X <- ( x@y[ , "11" ] + x@y[ , "10" ] ) == 1
         selection0X <- ( x@y[ , "01" ] + x@y[ , "00" ] ) == 1
         selectionX1 <- ( x@y[ , "11" ] + x@y[ , "01" ] ) == 1
         selectionX0 <- ( x@y[ , "10" ] + x@y[ , "00" ] ) == 1
      }

      # only considering the first probit equation
      result$IMR1X <- NA
      result$IMR1X[ selection1X ] <- dnorm( x@predictors[ selection1X, 1 ] ) /
         pnorm( x@predictors[ selection1X, 1 ] )

      result$delta1X <- NA
      result$delta1X[ selection1X ] <- result$IMR1X[ selection1X ] *
         ( result$IMR1X[ selection1X ] + x@predictors[ selection1X, 1 ] )

      result$IMR0X <- NA
      result$IMR0X[ selection0X ] <- -dnorm( x@predictors[ selection0X, 1 ] ) /
         pnorm( -x@predictors[ selection0X, 1 ] )
      result$delta0X <- NA
      result$delta0X[ selection0X ] <- result$IMR0X[ selection0X ] *
         (result$IMR0X[ selection0X ] + x@predictors[ selection0X, 1 ] )

      # only considering the second probit equation
      result$IMRX1 <- NA
      result$IMRX1[ selectionX1 ] <- dnorm( x@predictors[ selectionX1, 2 ] ) /
         pnorm( x@predictors[ selectionX1, 2 ] )

      result$deltaX1 <- NA
      result$deltaX1[ selectionX1 ] <- result$IMRX1[ selectionX1 ] *
         ( result$IMRX1[ selectionX1 ] + x@predictors[ selectionX1, 2 ] )

      result$IMRX0 <- NA
      result$IMRX0[ selectionX0 ] <- -dnorm( x@predictors[ selectionX0, 2 ] ) /
         pnorm( -x@predictors[ selectionX0, 2 ] )

      result$deltaX0 <- NA
      result$deltaX0[ selectionX0 ] <- result$IMRX0[ selectionX0 ] *
         ( result$IMRX0[ selectionX0 ] + x@predictors[ selectionX0, 2 ] )

   } else {
      stop( errorMessage )
   }
   result$no <- NULL
   return( result )
}
