cesInterN3 <- function( funcName, par, xNames, tName, data, rhoApprox ) {

      # interpolation if rho and/or rho_1 are close to zero
      coefArray <- array( par, c( length( par ), 2, 2 ) )
      dimnames( coefArray ) <- list( names( par ), 
         c( "rho_1 = 0", "rho_1 = E" ),  c( "rho = 0", "rho = E" ) )
      weights <- c( 0, 0 )
      names( weights ) <- c( "rho_1 = 0", "rho = 0" )
      rhoNames <- c( "rho_1", "rho" )
      for( i in 1:2 ) {
         if( abs( par[ rhoNames[ i ] ] ) <= rhoApprox ) {
            # permute the array so that the second dimension is for this 'i'
            atemp <- aperm( coefArray, c( 1, i + 1, 4 - i ) )
            atemp[ rhoNames[ i ], 1, ] <- 0
            atemp[ rhoNames[ i ], 2, ] <- 
               rhoApprox * (-1)^( par[ rhoNames[ i ] ] < 0 )
            # permute the array back to its initial state
            coefArray <- aperm( atemp, c( 1, i + 1, 4 - i ) )
            weights[ i ] <- 1 - abs( par[ rhoNames[ i ] ] ) / rhoApprox
         }
      }
      result <- 0
      weightMatrix <- cbind( weights, 1 - weights )
      for( i in 1:2 ) {
         for( j in 1:2 ) {
            if( weightMatrix[ 1, i ] != 0 && weightMatrix[ 2, j ] != 0 ) {
               result <- result + weightMatrix[ 1, i ] * weightMatrix[ 2, j ] * 
                  do.call( funcName, args = list( coef = coefArray[ , i, j ], 
                     data = data, xNames = xNames, tName = tName ) )
            }
         }
      }

   return( result )
}
