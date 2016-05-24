###This file contains useful functions not to be called by the user

## Initialize f
'initialf' <- 
  function(x,n,Xin,h){
    mean(apply(dnorm(rep(1,n)%o%x-Xin,sd = rep(1,n)%o%h),1,prod))
  }

## Find initial y
'initialy' <-
  function(Xin){
    n <- nrow(Xin)
    d <- ncol(Xin)
    h <- apply(Xin,2,sd)*n^(-1/(d+4))
    log(apply(Xin,1,initialf,n=n,Xin=Xin,h=h))
}

## Various functions for computing integral over unit simplex of
## exp( y[1]w[1] + y[2]w[2] + ... + y[d+1]( 1 - w[1] + ... + w[d] )

## Assuming y is an ordered vector
'JAD.ord' <-  function( y, d=length(y) - 1, eps=10^-3 ) {
  d <- length( y ) - 1
  if( y[ d + 1 ] - y[ 1 ] < eps  )
    return( JAD.approx( y, d ) )
  else {
    if( d==1 ) {
      return( ( exp( y[ 2 ] ) - exp( y[ 1 ] ) ) / (y[ 2 ] - y[ 1 ]) )
    }
    else {
      return( ( JAD.ord( y[ - 1], d-1, eps ) - JAD.ord( y[ - ( d + 1 ) ], d-1, eps ) ) / ( y [ d + 1 ]  - y[1] ) )
    }
  }
}

'JAD' <- function( y, eps=10^-3 ) {
  d <- length( y ) - 1
  return( JAD.ord( sort.int( y ) ,d ) )
}

'JAD.approx' <- function( y, d = length(y) - 1 ) {
  z <- y - mean( y )
  ans <- 1 + sum( z^2 ) / (2 * ( d + 1 ) * ( d + 2 ) )  + sum( z ^ 3 ) / ( 3 * ( d + 1 ) * ( d + 2 ) * ( d + 3 ) )
  ans <- ans / factorial( d ) * exp( mean( y ) )
  return( ans )
}

## Return the loglikelihood given proportions in mixture
'logliklcdmix' <- function( y, props ) {
  if( length( props ) == 1 ) {
    return( sum( y ) )
  }
  else {
    return( sum( log( apply( apply( exp( y ), 1, "*", props ), 2, sum)  ) ) )
  }
}
