## --- derivatives of the Hessian with respect to the betas ---
snqProfitHessianDeriv <- function( prices, weights, nFix = 0, form = 0 ) {
   prices   <- array( prices )
   weights <- array( weights )
   nNetput  <- dim(array(prices))
   nCoef    <- nNetput + nNetput * ( nNetput - 1 )/2 + nNetput * nFix
      # number of coefficients
   if( form == 0 ) {
      nCoef <- nCoef + ( nFix + 1 ) * nFix / 2
   } else if ( form == 1 ) {
      nCoef <- nCoef + nNetput * ( nFix + 1 ) * nFix / 2
   } else {
      stop( "argument 'form' must be either 0 or 1" )
   }
   normPrice <- sum( t( prices ) %*% weights )
   Hderiv    <- array( 0, c( nNetput * ( nNetput - 1 )/2, nCoef ) )
   kro       <- diag( 1, nNetput, nNetput )    #identity matrix for Kronecker Delta
   for( i in 1:( nNetput - 1 ) ) {
      for( j in i:( nNetput - 1 ) ) {
         for( k in 1:( nNetput - 1 ) ) {
            for( l in k:( nNetput - 1 ) ) {
               Hderiv[ veclipos( i, j, nNetput-1 ),
                  nNetput + veclipos( k, l, nNetput - 1 ) ] <-
                     ( kro[ i, k ] * kro[ j, l ] +
                     kro[ i, l ] * kro[ j, k ] * ( 1 - kro[ i, j ] ) ) /
                     normPrice -
                     kro[ i, l ] * weights[ j ] *
                     ( prices[ k ] - prices[ nNetput ] ) / normPrice^2 -
                     ( 1 - kro[ k, l ] ) * kro[ i, k ] * weights[ j ] *
                     ( prices[ l ] - prices[ nNetput ] ) / normPrice^2 -
                     kro[ j, k ] * weights[ i ] *
                     ( prices[ l ] - prices[ nNetput ] ) / normPrice^2 -
                     ( 1 - kro[ k, l ] ) * kro[ j, l ] * weights[ i ] *
                     ( prices[ k ] - prices[ nNetput ] ) / normPrice^2 +
                     ( 2 - kro[ k, l ] ) * weights[ i ] * weights[ j ] *
                     ( prices[ k ] - prices[ nNetput ] ) *
                     ( prices[ l ] - prices[ nNetput ] ) / normPrice^3
            }
         }
      }
   }
   return( Hderiv )
}
