cesEstKmenta <- function( yName, xNames, data, vrs ){

   result <- list()

   ## number of exogenous variables
   nExog <- length( xNames )

   ## Estimating the unrestricted translog model
   result$translog <- translogEst( yName = yName, xNames = xNames, data = data )

   ## Testing restrictions implied by the Kmenta approximation - Wald test
   restrictions <- c( "b_1_2 = -b_1_1", "b_1_2 = -b_2_2" )
   if( !vrs ) {
      restrictions <- c( "a_1 + a_2 = 1", restrictions )
   }
   result$testKmenta <- lht( result$translog$est, restrictions )

   ## Estimating restricted model
   result$kmenta <- systemfit( formula = formula(result$translog$est),
      data = model.frame( result$translog$est ),
      restrict.matrix = gsub( "([ab]\\_)", "eq1_\\1", restrictions) )

   ## Parameter vector
   result$coefficients <- numeric( 4 )
   names( result$coefficients ) <- cesCoefNames( nExog, vrs )

   ## Defining gamma
   result$coefficients[ "gamma" ] <- exp( coef( result$kmenta )[ "eq1_(Intercept)" ] )

   ## Defining nu
   result$coefficients[ "nu" ] <- coef( result$kmenta )[ "eq1_a_1" ] +
                             coef( result$kmenta )[ "eq1_a_2" ]

   ## Defining delta
   result$coefficients[ "delta" ] <- coef( result$kmenta )[ "eq1_a_1" ] /
                               result$coefficients[ "nu" ]

   ## Defining rho
   result$coefficients[ "rho" ] <-
      coef( result$kmenta )[ "eq1_b_1_2" ] * result$coefficients[ "nu" ] /
      ( coef( result$kmenta )[ "eq1_a_1" ] * coef( result$kmenta )[ "eq1_a_2" ] )

   ## Delta method
   jacobian <- matrix( 0, nrow = length( result$coefficients ),
      ncol = length( coef( result$kmenta )))
   rownames( jacobian ) <- names( result$coefficients )
   colnames( jacobian ) <- names( coef( result$kmenta ))

   jacobian[ "gamma", "eq1_(Intercept)" ] <-
      exp( coef( result$kmenta )[ "eq1_(Intercept)" ] )
   jacobian[ "nu", c( "eq1_a_1", "eq1_a_2" ) ] <- 1
   jacobian[ "delta", "eq1_a_1" ] <- 1 / result$coefficients[ "nu" ] -
      coef( result$kmenta )[ "eq1_a_1" ] / result$coefficients[ "nu" ]^2
   jacobian[ "delta", "eq1_a_2" ] <-
      - coef( result$kmenta )[ "eq1_a_1" ] / result$coefficients[ "nu" ]^2
   jacobian[ "rho", "eq1_a_1" ] <- coef( result$kmenta )[ "eq1_b_1_2" ] /
      ( coef( result$kmenta )[ "eq1_a_1" ] * coef( result$kmenta )[ "eq1_a_2" ] ) -
      coef( result$kmenta )[ "eq1_b_1_2" ] * result$coefficients[ "nu" ] /
      coef( result$kmenta )[ "eq1_a_1" ]^2 / coef( result$kmenta )[ "eq1_a_2" ]
   jacobian[ "rho", "eq1_a_2" ] <- coef( result$kmenta )[ "eq1_b_1_2" ] /
      ( coef( result$kmenta )[ "eq1_a_1" ] * coef( result$kmenta )[ "eq1_a_2" ] ) -
      coef( result$kmenta )[ "eq1_b_1_2" ] * result$coefficients[ "nu" ] /
      coef( result$kmenta )[ "eq1_a_1" ] / coef( result$kmenta )[ "eq1_a_2" ]^2
   jacobian[ "rho", "eq1_b_1_2" ] <-  result$coefficients[ "nu" ] /
      coef( result$kmenta )[ "eq1_a_1" ] / coef( result$kmenta )[ "eq1_a_2" ]
   result$vcov <- jacobian %*% vcov( result$kmenta ) %*%  t( jacobian )

   if( !vrs ) {
      result$coefficients <- result$coefficients[ 1:3 ]
      result$vcov <- result$vcov[ 1:3, 1:3 ]
   }

   return(result)
}

