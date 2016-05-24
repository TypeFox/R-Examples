# calculate part "B1"
cesDerivCoefN4B1 <- function( coef, data, xNames ) {

      B1 <- coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) + 
         ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ])

   return( B1 )
}


# calculate part "L1"
cesDerivCoefN4L1 <- function( coef, data, xNames ) {

   L1 <- coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) + 
      ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )

   return( L1 )
}


# calculate part "B2"
cesDerivCoefN4B2 <- function( coef, data, xNames ) {

      B2 <- coef[ "delta_2" ] * data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) + 
         ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ])

   return( B2 )
}


# calculate part "L2"
cesDerivCoefN4L2 <- function( coef, data, xNames ) {

   L2 <- coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) + 
      ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] )

   return( L2 )
}


# calculate part "B"
cesDerivCoefN4B <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

      B <- coef[ "delta" ] * B1^( coef[ "rho" ] / coef[ "rho_1" ] ) + 
         ( 1 - coef[ "delta" ] ) * B2^( coef[ "rho" ] / coef[ "rho_2" ] )

   return( B )
}


# derivatives with respect to gamma
cesDerivCoefN4Gamma <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN4L1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   L2 <- cesDerivCoefN4L2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- exp( - coef[ "nu" ] * ( - coef[ "delta" ] * L1 -
                  ( 1 - coef[ "delta" ] ) * L2 ) )
         } else {
            result <- exp( - coef[ "nu" ] * ( - coef[ "delta" ] * L1 +
                  ( 1- coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- exp( - coef[ "nu" ] * ( coef[ "delta" ] *
               log( B1 ) / coef[ "rho_1" ] -
               ( 1 - coef[ "delta" ] ) * L2 ) )
      } else {
         result <- 
            exp( - coef[ "nu" ] *
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- exp( - ( coef[ "nu" ] / coef[ "rho" ] ) *
            log( coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) +
            ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 )
            ) )
      } else {
         result <- 
            ( coef[ "delta" ] * 
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <-
         ( coef[ "delta" ] * B1^( coef[ "rho" ] / coef[ "rho_1" ] ) +
            ( 1 - coef[ "delta" ] ) * 
               exp( - coef[ "rho" ] * L2 ) 
         )^( - ( coef[ "nu" ] / coef[ "rho" ] ) )
   } else {
      result <- B^(-coef[ "nu" ]/coef[ "rho" ])
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}


# derivatives with respect to lambda
cesDerivCoefN4Lambda <- function( coef, data, xNames, tName ) {

   if( is.null( tName ) || ! "lambda" %in% names( coef ) ) {
      stop( "internal error: cannot calculate derivative w.r.t. lambda",
         " if 'tName' or 'lambda' is not given" )
   }
   
   result <- cesDerivCoefN4Gamma( coef = coef, data = data, 
      xNames = xNames, tName = tName ) * coef[ "gamma" ] * data[[ tName ]]

   return( result )
}


# derivatives with respect to delta_1
cesDerivCoefN4Delta1 <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN4L1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   L2 <- cesDerivCoefN4L2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] *
                  ( - log( data[[ xNames[ 1 ] ]] ) +
                     log( data[[ xNames[ 2 ] ]] ) ) *
               exp( -coef[ "nu" ] * 
                  ( - coef[ "delta" ] * L1 -
                     ( 1 - coef[ "delta" ] ) * L2 ) )
         } else {
            result <- coef[ "gamma" ] *
               ( - coef[ "nu" ]* 
               ( coef[ "delta" ] *
                  ( - log( data[[ xNames[ 1 ] ]] ) +
                     log( data[[ xNames[ 2 ] ]] ) ) ) ) *
               exp( -coef[ "nu" ] * 
                  ( - coef[ "delta" ] * L1 +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] * coef[ "delta" ] * 
               ( data[[ xNames[ 1 ] ]]^( - coef[ "rho_1" ] ) -
                  data[[ xNames[ 2 ] ]]^( - coef[ "rho_1" ] ) ) ) /
            ( coef[ "rho_1" ] * B1 ) *
            exp( -coef[ "nu" ] *
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] -
                  ( 1 - coef[ "delta" ] ) * L2 ) )
      } else {
         result <- coef[ "gamma" ] *
            ( - coef[ "nu" ] / coef[ "rho_1" ] ) * 
            ( coef[ "delta" ] *
               ( data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
                  data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) / B1 ) *
            exp( -coef[ "nu" ] * 
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 )
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) *
                  coef[ "rho" ] * ( - log( data[[ xNames[ 1 ] ]] ) + 
                     log( data[[ xNames[ 2 ] ]] ) ) )
      } else {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) *
            coef[ "rho" ] *
            ( - log( data[[ xNames[ 1 ] ]] ) + log( data[[ xNames[ 2 ] ]] ) ) )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( ( 1 - coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( coef[ "delta" ] * coef[ "rho" ] *
            B1^( coef[ "rho" ] / coef[ "rho_1" ] - 1 ) *
            ( data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
               data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) /
            coef[ "rho_1" ] )
   } else {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
         (coef[ "rho" ]/coef[ "rho_1" ]) * coef[ "delta" ] * 
         B1^((coef[ "rho" ]-coef[ "rho_1" ])/coef[ "rho_1" ]) * 
         ( data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
            data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}


# derivatives with respect to delta_2
cesDerivCoefN4Delta2 <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN4L1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   L2 <- cesDerivCoefN4L2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- - coef[ "gamma" ] * coef[ "nu" ] * 
               ( 1 - coef[ "delta" ] ) *
               ( - log( data[[ xNames[ 3 ] ]] ) +
                  log( data[[ xNames[ 4 ] ]] ) ) *
               exp( -coef[ "nu" ] * 
                  ( - coef[ "delta" ] * L1 -
                     ( 1 - coef[ "delta" ] ) * L2 ) )
         } else {
            result <- coef[ "gamma" ] * 
               ( -coef[ "nu" ] * ( 1 - coef[ "delta" ] ) * 
                  ( data[[ xNames[ 3 ] ]]^( - coef[ "rho_2" ] ) -
                     data[[ xNames[ 4 ] ]]^( - coef[ "rho_2" ] ) ) ) /
               ( coef[ "rho_2" ] * B2 ) *
               exp( -coef[ "nu" ] *
                  ( - coef[ "delta" ] * L1 +
                     ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] *
            ( - coef[ "nu" ]* 
            ( ( 1 - coef[ "delta" ] ) *
               ( - log( data[[ xNames[ 3 ] ]] ) +
                  log( data[[ xNames[ 4 ] ]] ) ) ) ) *
            exp( -coef[ "nu" ] * 
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] -
               ( 1 - coef[ "delta" ] ) * L2 ) )
      } else {
         result <- coef[ "gamma" ] *
            ( - coef[ "nu" ] / coef[ "rho_2" ] ) * 
            ( ( 1 - coef[ "delta" ] ) *
               ( data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) - 
                  data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ]) ) / B2 ) *
            exp( -coef[ "nu" ] * 
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- - coef[ "gamma" ] * ( coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 )
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 ) *
                  coef[ "rho" ] * ( - log( data[[ xNames[ 3 ] ]] ) + 
                     log( data[[ xNames[ 4 ] ]] ) ) )
      } else {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( ( 1 - coef[ "delta" ] ) * coef[ "rho" ] *
               B2^( coef[ "rho" ] / coef[ "rho_2" ] - 1 ) *
               ( data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) - 
                  data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ]) ) /
               coef[ "rho_2" ] )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( ( 1 - coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( ( 1 - coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) *
         coef[ "rho" ] *
         ( - log( data[[ xNames[ 3 ] ]] ) + log( data[[ xNames[ 4 ] ]] ) ) )
   } else {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
         (coef[ "rho" ]/coef[ "rho_2" ]) * ( 1 - coef[ "delta" ] ) * 
         B2^((coef[ "rho" ]-coef[ "rho_2" ])/coef[ "rho_2" ]) * 
         ( data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) - 
            data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ]) )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}

# derivatives with respect to delta
cesDerivCoefN4Delta <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN4L1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   L2 <- cesDerivCoefN4L2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- - coef[ "gamma" ] * coef[ "nu" ] *
               ( - L1 + L2 ) * 
               exp( - coef[ "nu" ] * 
                  ( - coef[ "delta" ] * L1 -
                     ( 1 - coef[ "delta" ] ) * L2 ) )
         } else {
            result <- - coef[ "gamma" ] * coef[ "nu" ] *
               ( - L1 - log( B2 ) / coef[ "rho_2" ] ) *
               exp( - coef[ "nu" ] * 
                  ( - coef[ "delta" ] * L1 +
                     ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- - coef[ "gamma" ] * coef[ "nu" ] *
            ( log( B1 ) / coef[ "rho_1" ] + L2 ) *
            exp( - coef[ "nu" ] * 
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] -
                  ( 1 - coef[ "delta" ] ) * L2 ) )
      } else {
         result <- coef[ "gamma" ] *
            ( - coef[ "nu" ] * 
               ( log( B1 ) / coef[ "rho_1" ] - log( B2 ) / coef[ "rho_2" ] ) ) *
            exp( -coef[ "nu" ] * 
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( exp( - coef[ "rho" ] * L1 ) -
               exp( - coef[ "rho" ] * L2 )
            ) 
      } else {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( exp( - coef[ "rho" ] * L1 ) -
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( ( 1 -coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( - exp( - coef[ "rho" ] * L2 ) +
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) )
   } else {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
         ( B1^(coef[ "rho" ]/coef[ "rho_1" ]) - 
            B2^(coef[ "rho" ]/coef[ "rho_2" ]) )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}


# derivatives with respect to rho_1
cesDerivCoefN4Rho1 <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN4L1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   L2 <- cesDerivCoefN4L2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] * 
               ( 0.5 * 
                  ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] )^2 +
                     ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )^2 ) -
                  0.5 * L1^2 ) *
               exp( - coef[ "nu" ] *
                  ( - coef[ "delta" ] * L1 -
                     ( 1 - coef[ "delta" ] ) * L2 ) )
         } else {
            result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] * 
               ( 0.5 * 
                  ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] )^2 +
                     ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )^2 ) -
                  0.5 * L1^2 ) *
               exp( - coef[ "nu" ] *
                  ( - coef[ "delta" ] * L1 +
                     ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] * 
            ( - log( B1 ) / coef[ "rho_1" ]^2 +
               ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
                  data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ] ) - 
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) * 
                  data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ] ) ) /
                  ( coef[ "rho_1" ] * B1 ) ) *
            exp( - coef[ "nu" ] *
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] -
                  ( 1 - coef[ "delta" ] ) * L2 ) )
      } else {
         result <- coef[ "gamma" ] * 
            ( - coef[ "nu" ] *
               ( coef[ "delta" ] * coef[ "rho_1" ] *
                  ( -coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
                     data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ] ) - 
                     ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) * 
                     data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ] ) ) /
                     ( B1 * coef[ "rho_1" ]^2   ) -
                     coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ]^2 ) ) *
            exp( - coef[ "nu" ] *
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- - coef[ "gamma" ] * coef[ "delta" ] * coef[ "nu" ] * 
            ( coef[ "delta" ] *  exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 )
            )^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            exp( - coef[ "rho" ] * L1 ) *
                  ( - 0.5 * L1^2 +
                     0.5 * ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] )^2 +
                        ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )^2 )
                  )
      } else {
         result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] * 
            ( coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] )
            )^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            exp( - coef[ "rho" ] * L1 ) *
            ( - 0.5 * L1^2 +
               0.5 * ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] )^2 +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )^2 )
            )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         ( ( 1 - coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] )
         )^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( coef[ "delta" ] * log( B1 ) * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) *
            ( - coef[ "rho" ] / coef[ "rho_1" ]^2 ) +
            coef[ "delta" ] * ( coef[ "rho" ] / coef[ "rho_1" ] ) *
            B1^( coef[ "rho" ] / coef[ "rho_1" ] - 1 ) *
            ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) *
               data[[ xNames[ 1 ] ]]^( - coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) *
               data[[ xNames[ 2 ] ]]^( - coef[ "rho_1" ] )
            )
         )
   } else {
         result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) * 
            B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
            ( coef[ "delta" ] * log( B1 ) * B1^(coef[ "rho" ]/coef[ "rho_1" ]) * 
               ( -coef[ "rho" ]/coef[ "rho_1" ]^2 ) + 
               coef[ "delta" ] * 
               B1^((coef[ "rho" ]-coef[ "rho_1" ])/coef[ "rho_1" ]) * 
               (coef[ "rho" ]/coef[ "rho_1" ]) *
               ( -coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
                  data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) * 
                  data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}


# derivatives with respect to rho_2
cesDerivCoefN4Rho2 <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN4L1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   L2 <- cesDerivCoefN4L2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- - coef[ "gamma" ] * coef[ "nu" ] * ( 1 - coef[ "delta" ] ) * 
               ( 0.5 * 
                  ( coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] )^2 +
                     ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] )^2 ) -
                  0.5 * L2^2 ) *
               exp( - coef[ "nu" ] *
                  ( - coef[ "delta" ] * L1 -
                     ( 1 - coef[ "delta" ] ) * L2 ) )
         } else {
            result <- - coef[ "gamma" ] * coef[ "nu" ] *
               ( 1 - coef[ "delta" ] ) *
               ( - log( B2 ) / coef[ "rho_2" ]^2 +
                  ( -coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) * 
                     data[[ xNames[ 3 ] ]]^( - coef[ "rho_2" ]) - 
                     ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) * 
                     data[[ xNames[ 4 ] ]]^( - coef[ "rho_2" ] ) ) /
                  ( coef[ "rho_2" ] * B2 ) ) *
               exp( - coef[ "nu" ] *
                  ( - coef[ "delta" ] * L1 +
                     ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- - coef[ "gamma" ] * coef[ "nu" ] * ( 1 - coef[ "delta" ] ) *
            ( 0.5 * 
               ( coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] )^2 +
                  ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] )^2 ) -
               0.5 * L2^2 ) *
            exp( - coef[ "nu" ] *
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] -
                  ( 1 - coef[ "delta" ] ) * L2 ) )
      } else {
         result <- coef[ "gamma" ] * 
            ( - coef[ "nu" ] *
               ( ( 1 - coef[ "delta" ] ) * coef[ "rho_2" ] *
                  ( -coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) * 
                     data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ] ) - 
                     ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) * 
                     data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ] ) ) /
                     ( B2 * coef[ "rho_2" ]^2   ) -
                     ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ]^2 ) ) *
            exp( - coef[ "nu" ] *
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- - coef[ "gamma" ] * ( 1 - coef[ "delta" ] ) * coef[ "nu" ] *
            ( ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 ) +
               coef[ "delta" ] * exp( - coef[ "rho" ] * L1 )
            )^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            exp( - coef[ "rho" ] * L2 ) *
               ( - 0.5 * L2^2 +
                  0.5 * ( coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] )^2 +
                     ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] )^2 )
               )
      } else {
         result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( ( 1 - coef[ "delta" ] ) * log( B2 ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) *
               ( - coef[ "rho" ] / coef[ "rho_2" ]^2 ) +
               ( 1 - coef[ "delta" ] ) * ( coef[ "rho" ] / coef[ "rho_2" ] ) *
               B2^( coef[ "rho" ] / coef[ "rho_2" ] - 1 ) *
               ( -coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) * 
                  data[[ xNames[ 3 ] ]]^( - coef[ "rho_2" ]) - 
                  ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) * 
                  data[[ xNames[ 4 ] ]]^( - coef[ "rho_2" ] ) ) )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- - coef[ "gamma" ] * coef[ "nu" ] * ( 1 - coef[ "delta" ] ) * 
         ( ( 1 - coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] )
         )^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         exp( - coef[ "rho" ] * L2 ) *
         ( - 0.5 * L2^2 +
            0.5 * ( coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] )^2 +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] )^2 )
         )
   } else {
         result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) * 
            B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
            ( ( 1 - coef[ "delta" ] ) * log( B2 ) * 
               B2^(coef[ "rho" ]/coef[ "rho_2" ]) * 
               ( -coef[ "rho" ]/coef[ "rho_2" ]^2 ) + 
               ( 1 - coef[ "delta" ] ) * 
               B2^((coef[ "rho" ]-coef[ "rho_2" ])/coef[ "rho_2" ]) * 
               (coef[ "rho" ]/coef[ "rho_2" ]) *
               ( -coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) * 
                  data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) - 
                  ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) * 
                  data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ]) ) )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}

# derivatives with respect to rho
cesDerivCoefN4Rho <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN4L1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   L2 <- cesDerivCoefN4L2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- coef[ "gamma" ] * coef[ "nu" ] *
               ( -0.5 * ( coef[ "delta" ] * L1^2 +
                  ( 1 - coef[ "delta" ] ) * L2^2 ) +
                  0.5 * ( coef[ "delta" ] * L1 +
                     ( 1 - coef[ "delta" ] ) * L2 )^2 ) *
               exp( - coef[ "nu" ] * 
                  ( - coef[ "delta" ] * L1 -
                     ( 1 - coef[ "delta" ] ) * L2 ) )
         } else {
            result <- coef[ "gamma" ] * coef[ "nu" ] *
               ( -0.5 * ( coef[ "delta" ] * L1^2 +
                  ( 1 - coef[ "delta" ] ) * ( log( B2 ) / coef[ "rho_2" ] )^2 ) +
                  0.5 * ( - coef[ "delta" ] * L1 +
                     ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] )^2 ) *
               exp( - coef[ "nu" ] * 
                  ( - coef[ "delta" ] * L1 +
                     ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else  if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] * coef[ "nu" ] *
            ( -0.5 * ( coef[ "delta" ] * ( log( B1 ) / coef[ "rho_1" ] )^2 +
               ( 1 - coef[ "delta" ] ) * L2^2 ) +
               0.5 * ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] -
                  ( 1 - coef[ "delta" ] ) * L2 )^2 ) *
            exp( - coef[ "nu" ] * 
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] -
                  ( 1 - coef[ "delta" ] ) * L2 ) )
      } else {
         result <- coef[ "gamma" ] * coef[ "nu" ] *
            ( -0.5 * ( coef[ "delta" ] * log( B1 )^2 / coef[ "rho_1" ]^2 +
                  ( 1 - coef[ "delta" ] ) * log( B2 )^2 / coef[ "rho_2" ]^2 ) +
               0.5 * ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] )^2 ) *
            exp( - coef[ "nu" ] * 
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] *
            log( coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) + 
               ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 )
            ) *
            ( coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "nu" ] / coef[ "rho" ]^2 ) -
            coef[ "gamma" ] * ( coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) * 
            ( - coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) * L1 -
               ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 ) *
               L2 )
      } else {
         result <- coef[ "gamma" ] *
            log( coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) ) *
            ( coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "nu" ] / coef[ "rho" ]^2 ) +
            coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] )
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( - coef[ "delta" ] *
               exp( - coef[ "rho" ] * L1 ) * L1 +
               ( 1 - coef[ "delta" ] ) * log( B2 ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) / coef[ "rho_2" ] )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] *
         log( ( 1 - coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] *
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) ) *
         ( ( 1 - coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] ) *
         ( coef[ "nu" ] / coef[ "rho" ]^2 ) +
         coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( ( 1 - coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] *
               B1^( coef[ "rho" ] / coef[ "rho_1" ] )
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( - ( 1 - coef[ "delta" ] ) *
            exp( - coef[ "rho" ] * L2 ) * L2 +
            ( coef[ "delta" ] * log( B1 ) * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) / coef[ "rho_1" ] )
         )
   } else {
         result <- coef[ "gamma" ] * log( B ) * 
            B^(-coef[ "nu" ]/coef[ "rho" ]) * ( coef[ "nu" ] / coef[ "rho" ]^2 ) +
            coef[ "gamma" ] * ( -coef[ "nu" ]/coef[ "rho" ] ) * 
            B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) * 
            ( coef[ "delta" ] * log( B1 ) * 
               B1^(coef[ "rho" ]/coef[ "rho_1" ]) / coef[ "rho_1" ] + 
               ( 1 - coef[ "delta" ] ) * log( B2 ) * 
               B2^(coef[ "rho" ]/coef[ "rho_2" ]) / coef[ "rho_2" ] )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}

# derivatives with respect to nu
cesDerivCoefN4Nu <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN4L1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   L2 <- cesDerivCoefN4L2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- coef[ "gamma" ] *
               ( coef[ "delta" ] * L1 +
                  ( 1 - coef[ "delta" ] ) * L2 ) *
               exp( - coef[ "nu" ] *
                  ( - coef[ "delta" ] * L1 -
                     ( 1 - coef[ "delta" ] ) * L2 ) )
         } else {
            result <- coef[ "gamma" ] *
               ( coef[ "delta" ] * L1 -
               ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) *
               exp( - coef[ "nu" ] *  
                  ( - coef[ "delta" ] * L1 +
                     ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] *
            ( - coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta" ] ) * L2 ) *
            exp( - coef[ "nu" ] *  
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] -
                  ( 1 - coef[ "delta" ] ) * L2 ) )
      } else {
         result <- - coef[ "gamma" ] * 
            ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) *
            exp( - coef[ "nu" ] *
               ( coef[ "delta" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- - ( coef[ "gamma" ] / coef[ "rho" ] ) * 
            log( ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 ) +
               coef[ "delta" ] * exp( - coef[ "rho" ] * L1 )
            ) *            
            ( coef[ "delta" ] * exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * exp( - coef[ "rho" ] * L2 )
            )^( - coef[ "nu" ] / coef[ "rho" ] )
      } else {
         result <- - ( coef[ "gamma" ] / coef[ "rho" ] ) * 
            log( coef[ "delta" ] * 
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) ) *
            ( coef[ "delta" ] * 
               exp( - coef[ "rho" ] * L1 ) +
               ( 1 - coef[ "delta" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- - ( coef[ "gamma" ] / coef[ "rho" ] ) * 
         log( ( 1 - coef[ "delta" ] ) * 
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) ) *
         ( ( 1 - coef[ "delta" ] ) * 
            exp( - coef[ "rho" ] * L2 ) +
            coef[ "delta" ] * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   } else {
         result <- - coef[ "gamma" ] * log( B ) * 
            B^( -coef[ "nu" ] / coef[ "rho" ] ) / coef[ "rho" ]
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}

