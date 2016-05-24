## ------- calculation of Hessian -------------
snqProfitHessian <- function( beta, prices, weights,
      scalingFactors = rep( 1, length( weights ) ) ) {

   prices <- unlist( prices ) * scalingFactors
   normPrice <- sum( t( prices ) %*% weights )
   hessian <- beta / normPrice -
      beta %*% prices %*% t( weights ) / normPrice^2 -
      weights %*% t( prices ) %*% beta / normPrice^2 +
      weights %*% t( weights ) *
      mean( ( t( prices ) %*% beta %*% prices ) / normPrice^3 )
   return( hessian )
}
