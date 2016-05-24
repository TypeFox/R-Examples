translogCheckMono <- function( xNames, data, coef, increasing = TRUE,
   strict = FALSE, dataLogged = FALSE,
   tol = 10 * .Machine$double.eps ) {

   if( !is.logical( increasing ) ) {
      stop( "argument 'increasing' must be logical" )
   } else if( length( increasing ) == 1 ) {
      increasing <- rep( increasing, length( xNames ) )
   } else if( length( increasing ) != length( xNames ) ) {
      stop( "argument 'increasing' must be either a single logical variable",
         " or a vector of logical variables with same length as argument",
         " 'xNames' (", length( xNames ), ")" )
   }

   result <- list()

   # strip coefficients of shifter variables
   coef <- coef[ !grepl( "^d_[0-9]*$", names( coef ) ) ]
   coef <- coef[ !grepl( "^d_[0-9]*_", names( coef ) ) ]

   deriv <- translogDeriv( xNames = xNames, data = data, coef = coef,
      dataLogged = dataLogged )$deriv

   nExog <- ncol( deriv )
   nObs <- nrow( deriv )

   result$exog <- matrix( NA, nrow = nObs, ncol = nExog )
   colnames( result$exog ) <- colnames( deriv )

   for( i in 1:nExog ) {
      if( strict ) {
         result$exog[ , i ] <- deriv[ , i ] * (-1)^increasing[ i ] < 0
      } else {
         result$exog[ , i ] <- deriv[ , i ] * (-1)^increasing[ i ] <= tol
      }
   }

   result$obs <- rowSums( !result$exog ) == 0
   result$increasing <- increasing
   result$strict     <- strict

   class( result ) <- "translogCheckMono"
   return( result )
}

