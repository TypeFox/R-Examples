cesCoefAddRho <- function( coef, vrs, rho1, rho2, rho, nExog, nested ) {

   if( !is.null( rho ) ) {
      if( vrs ) {
         coef <- c( coef[ -length( coef ) ], rho, coef[ length( coef ) ] )
      } else {
         coef <- c( coef, rho )
      }
      if( !is.null( names( coef ) ) ) {
         names( coef )[ length( coef ) - vrs ] <- "rho"
      }
   }
   if( !is.null( rho2 ) && nested && nExog == 4 ) {
      coefAfter <- ( length( coef ) - vrs ):( length( coef ) )
      coef <- c( coef[ -coefAfter ], rho2, coef[ coefAfter ] )
      if( !is.null( names( coef ) ) ) {
         names( coef )[ coefAfter[1] ] <- "rho_2"
      }
   }
   if( !is.null( rho1 ) && nested ) {
      coefAfter <- ( length( coef ) - vrs - ( nExog == 4 ) ):( length( coef ) )
      coef <- c( coef[ -coefAfter ], rho1, coef[ coefAfter ] )
      if( !is.null( names( coef ) ) ) {
         names( coef )[ coefAfter[1] ] <- "rho_1"
      }
   }
   return( coef )
}

