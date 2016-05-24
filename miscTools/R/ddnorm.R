## derivatives of the density function of the normal distribution
## with respect to x
ddnorm <- function( x, mean = 0, sd = 1 ) {
   deriv <- - dnorm( x = x, mean = mean, sd = sd ) * ( x - mean ) / sd^2
   return( deriv )
}
