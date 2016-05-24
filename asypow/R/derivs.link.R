derivs.link <- function(u, link){

### Value and first two derivatives of link function

### link
###      (1) Logistic
###          g(u) = exp(u)/(1+exp(u))
###      (2) Complementary log-log
###          g(u) = 1 - exp( -exp(u) )

### Returns vector of length three 
###     [1] - value of link function at u
###     [2] - first derivative of link function wrt u
###     [3] - second derivative of link function wrt u

  ans <- rep( 0, 3 )
  ans[1] <- 1
  BIG <- 690
  if ( u > BIG ) return( ans )
  if (link < 1 || link > 2) stop('Illegal value of link')

### Logistic

  if (link == 1 ) {
    tmp <- exp(u)
    ans[1] <- tmp / ( 1 + tmp )
    ans[2] <- ans[1] * ( 1 - ans[1] )
    ans[3] <- ans[2] * ( 1 - 2*ans[1] )  } else {
      tmp <- exp(u)
      ans[1] <- 1 - exp( -tmp )
      ans[2] <- exp( u - tmp )
      ans[3] <- ans[2] - exp( 2*u - tmp )   }

  return( ans )     }
