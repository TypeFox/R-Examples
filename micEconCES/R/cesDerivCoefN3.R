# calculate part "B1"
cesDerivCoefN3B1 <- function( coef, data, xNames ) {

   result <- coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
      ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] )

   return( result )
}


# calculate part "L1"
cesDerivCoefN3L1 <- function( coef, data, xNames ) {

   L1 <- coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) + 
      ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )

   return( L1 )
}


# calculate part "B"
cesDerivCoefN3B <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   result <- coef[ "delta" ] *
      B1^( coef[ "rho" ] / coef[ "rho_1" ] ) + 
      ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] )

   return( result )
}


# derivatives with respect to gamma
cesDerivCoefN3Gamma <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN3L1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- exp( coef[ "nu" ] * 
            ( coef[ "delta" ] * L1 +
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <-
            exp( coef[ "nu" ] * 
               ( coef[ "delta" ] * 
                  ( - log( B1 ) / coef[ "rho_1" ] ) +
                  ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <-
         ( coef[ "delta" ] *
            exp( L1 )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   } else {
      result <- B^( - coef[ "nu" ] / coef[ "rho" ] )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}


# derivatives with respect to lambda
cesDerivCoefN3Lambda <- function( coef, data, xNames, tName ) {

   if( is.null( tName ) || ! "lambda" %in% names( coef ) ) {
      stop( "internal error: cannot calculate derivative w.r.t. lambda",
         " if 'tName' or 'lambda' is not given" )
   }
   
   result <- cesDerivCoefN3Gamma( coef = coef, data = data, 
      xNames = xNames, tName = tName ) * coef[ "gamma" ] * data[[ tName ]]

   return( result )
}

# derivatives with respect to delta_1
cesDerivCoefN3Delta1 <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN3L1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] *
               ( - log( data[[ xNames[ 1 ] ]] ) + log( data[[ xNames[ 2 ] ]] ) ) *
               exp( - coef[ "nu" ] * 
                  ( - coef[ "delta" ] * L1 -
                     ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) 
                  ) )
      } else {
         result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] *
            ( data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
               data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) /
            ( coef[ "rho_1" ] * B1 ) *
            exp( - coef[ "nu" ] * ( - coef[ "delta" ] * 
                  ( - log( B1 ) / coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <-
         coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] * 
         exp( L1 )^( -coef[ "rho" ] ) *
         ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) ) *
         ( coef[ "delta" ] *
            exp( L1 )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 )
   } else {
      result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] * 
         B1^( coef[ "rho" ] / coef[ "rho_1" ] - 1 ) / coef[ "rho_1" ] * 
         B^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) 
            - data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}


# derivatives with respect to delta
cesDerivCoefN3Delta <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN3L1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- - coef[ "gamma" ] * coef[ "nu" ] * 
            ( - L1 +
               log( data[[ xNames[ 3 ] ]] ) ) *
            exp( coef[ "nu" ] * ( coef[ "delta" ] * L1 +
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- - coef[ "gamma" ] * coef[ "nu" ] *
            ( log( B1 ) / coef[ "rho_1" ] + 
               log( data[[ xNames[ 3 ] ]] ) ) *
            exp( coef[ "nu" ] * ( coef[ "delta" ] * 
                  ( - log( B1 ) / coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <-
         - coef[ "gamma" ] * ( coef[ "nu" ] / coef[ "rho" ] ) *
         ( 
            exp( L1 )^( -coef[ "rho" ] ) -
            data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) *
         ( coef[ "delta" ] *
            exp( L1 )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 )
   } else {
      result <- -coef[ "gamma" ] * ( coef[ "nu" ] / coef[ "rho" ] ) *
         (
         B1^( coef[ "rho" ] / coef[ "rho_1" ] ) - 
         data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) *
         B^( -coef[ "nu" ] / coef[ "rho" ] - 1 )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}


# derivatives with respect to rho_1
cesDerivCoefN3Rho1 <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN3L1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] * 
            0.5 * ( coef[ "delta_1" ] * ( log( data[[ xNames[ 1 ] ]] ) )^2 +
               ( 1 - coef[ "delta_1" ] ) * ( log( data[[ xNames[ 2 ] ]] ) )^2 -
               L1^2 ) *
            exp( coef[ "nu" ] * ( coef[ "delta" ] * L1 +
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- - coef[ "gamma" ] * coef[ "nu" ] * 
            ( coef[ "delta" ] / coef[ "rho_1" ] ) *
            ( - log( B1 ) / coef[ "rho_1" ] + 
               ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
                  data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) -
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) *
                  data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) / B1 ) *
            exp( - coef[ "nu" ] * ( - coef[ "delta" ] * 
                  ( - log( B1 ) / coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] * 
         exp( -coef[ "rho" ] * L1 ) *
         0.5 * ( ( coef[ "delta_1" ] * ( log( data[[ xNames[ 1 ] ]] ) )^2 +
               ( 1 - coef[ "delta_1" ] ) * ( log( data[[ xNames[ 2 ] ]] ) )^2 ) -
            L1^2 ) *
         ( coef[ "delta" ] *
            exp( -coef[ "rho" ] * L1 ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( -coef[ "nu" ]/coef[ "rho" ] - 1 )
   } else {
      result <- -coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta" ] / 
         coef[ "rho_1" ] *
         B1^( coef[ "rho" ] / coef[ "rho_1" ] ) * 
         B^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( - log( B1 ) / coef[ "rho_1" ] +
         ( -coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
            data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) -
            ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) *
            data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) ) / B1 )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}


# derivatives with respect to rho
cesDerivCoefN3Rho <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN3L1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- coef[ "gamma" ] * coef[ "nu" ] *
            ( -0.5 * 
               ( coef[ "delta" ] * L1^2 +
                  ( 1- coef[ "delta" ] ) * ( log( data[[ xNames[ 3 ] ]] ) )^2 
               ) +
               0.5 * 
               ( - coef[ "delta" ] * L1 -
                  ( 1- coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) )^2 
            ) *
            exp( coef[ "nu" ] * ( coef[ "delta" ] * L1 +
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- - coef[ "gamma" ] * coef[ "nu" ] *
            ( 0.5 * ( coef[ "delta" ] * 
                  ( log( B1 ) / coef[ "rho_1" ] )^2 +
                  ( 1 - coef[ "delta" ] ) * ( log( data[[ xNames[ 3 ] ]] ) )^2 ) -
               0.5 * ( coef[ "delta" ] * 
                  (  log( B1 ) / coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) )^2 ) *
            exp( coef[ "nu" ] * ( - coef[ "delta" ] * 
                  ( log( B1 ) / coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- coef[ "gamma" ] * ( 
         ( coef[ "nu" ] / coef[ "rho" ]^2 ) * 
         log( coef[ "delta" ] *
            ( exp( L1 ) 
            )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta" ] ) * 
            data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) *
         ( coef[ "delta" ] *
            ( exp( L1 ) )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] ) -
         ( coef[ "nu" ] / coef[ "rho" ] ) * 
         ( coef[ "delta" ] *
            ( exp( L1 ) )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( -
            coef[ "delta" ] * L1 *
            ( exp( L1 ) )^( -coef[ "rho" ] ) -
            ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) *
               data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) )
   } else {
      result <- -coef[ "gamma" ] * ( coef[ "nu" ] / coef[ "rho" ] ) * 
         B^( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( -log( B ) / coef[ "rho" ] +
            B^(-1) * ( coef[ "delta" ] *
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) *
               ( log( B1 ) / coef[ "rho_1" ] ) -
            ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) *
            data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) )
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}


# derivatives with respect to nu
cesDerivCoefN3Nu <- function( coef, data, xNames, tName ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   L1 <- cesDerivCoefN3L1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- coef[ "gamma" ] *
            ( coef[ "delta" ] * L1 +
               ( 1- coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) *
            exp( - coef[ "nu" ] * ( - coef[ "delta" ] * L1 -
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- coef[ "gamma" ] * 
            ( coef[ "delta" ] *
               ( - log( B1 ) / coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) *
            exp( - coef[ "nu" ] * ( - coef[ "delta" ] * 
                  ( - log( B1 ) / coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- - ( coef[ "gamma" ] / coef[ "rho" ] ) *
         ( coef[ "delta" ] *
            exp( L1 )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] )
         )^( - coef[ "nu" ] / coef[ "rho" ] ) *
         log( coef[ "delta" ] *
            exp( L1 )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) )
   } else {
      result <- - coef[ "gamma" ] * log( B ) * 
         B^( -coef[ "nu" ] / coef[ "rho" ] ) / coef[ "rho" ]
   }

   if( !is.null( tName ) ){
      result <- result * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   return( result )
}

