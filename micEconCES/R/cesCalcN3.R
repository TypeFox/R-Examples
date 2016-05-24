cesCalcN3 <- function( xNames, tName, data, coef ) {

   gammaStar <- coef[ "gamma" ]
   if( "lambda" %in% names( coef ) ) {
      gammaStar <- gammaStar * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- gammaStar *
            exp( coef[ "nu" ] * 
               ( coef[ "delta" ] * (
                  coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) +
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- gammaStar *
            exp( coef[ "nu" ] * ( coef[ "delta" ] *
               ( -
                  log( coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
                     ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) ) /
                  coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- gammaStar * 
         ( coef[ "delta" ] *
            exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )
               )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   } else {
      result <-
         gammaStar * (
            coef[ "delta" ] *
            ( coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] )
            )^( coef[ "rho" ] / coef[ "rho_1" ] ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] )
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   }

   return( result )
}
