### Test the 'finalHessian' argument of optimization routines

library(maxLik)
set.seed( 4 )

                           # log-likelihood function, gradient, and Hessian for 1-parameter case (exponential distribution)
ll1i <- function(theta) {
   if(!all(theta > 0))
       return(NA)
   log(theta) - theta*t
}
ll1 <- function(theta) sum( log(theta) - theta*t )
gr1i <- function(theta) 1/theta - t
gr1 <- function(theta) sum( 1/theta - t )
hs1 <- function(theta) -100/theta^2
t <- rexp( 100, 2 )

## the same functions for 2-variable case (normal distribution)
ll2 <- function( param ) {
   ## log likelihood function
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   return( llValue )
}

## log likelihood function (individual observations)
ll2i <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   return( llValues )
}

gr2 <- function( param ) {
   ## function to calculate analytical gradients
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llGrad <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   return( llGrad )
}

## function to calculate analytical gradients (individual observations)
gr2i <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   llGrads <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   return( llGrads )
}

## function to calculate analytical Hessians
hs2 <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llHess <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llHess )
}
x <- rnorm(100, 1, 2)


## NR
# Estimate with only function values (single parameter)
a <- maxLik( ll1i, gr1i, start = 1, method = "NR" )
summary(a )
b <- maxLik( ll1i, gr1i, start = 1, method = "NR", finalHessian="bhhh")
                           # should issue a warning as BHHH not possible
summary(b )
c <- maxLik( ll1i, gr1i, start = 1, method = "NR", finalHessian=FALSE)
summary(c)
## (vector parameter)
a <- maxLik( ll2, gr2, start = c(0,1), method = "NR" )
summary(a )
b <- maxLik( ll2, gr2, start = c(0,1), method = "NR", finalHessian="bhhh")
                           # should issue a warning as BHHH not possible
summary(b )
c <- maxLik( ll2, gr2, start = c(0,1), method = "NR", finalHessian=FALSE)
summary(c)

## BFGSR
# Estimate with only function values (single parameter)
a <- maxLik( ll1i, gr1i, start = 1, method = "BFGSR" )
summary(a )
b <- maxLik( ll1i, gr1i, start = 1, method = "BFGSR", finalHessian="bhhh")
                           # should issue a warning as BHHH not possible
summary(b )
c <- maxLik( ll1i, gr1i, start = 1, method = "BFGSR", finalHessian=FALSE)
summary(c)
# Estimate with only function values (vector parameter)
a <- maxLik( ll2, gr2, start = c(0,1), method = "BFGSR" )
summary(a )
b <- maxLik( ll2, gr2, start = c(0,1), method = "BFGSR", finalHessian="bhhh")
                           # should issue a warning as BHHH not possible
summary(b )
c <- maxLik( ll2, gr2, start = c(0,1), method = "BFGSR", finalHessian=FALSE)
summary(c)


### Nelder-Mead
## Individual observations only
b <- maxLik( ll2i, start = c(0,1), method = "NM", finalHessian="bhhh")
summary(b)
## Individual observations, summed gradient
b <- maxLik( ll2i, gr2, start = c(0,1), method = "NM", finalHessian="bhhh")
                           # should issue a warning as BHHH not selected
                           # (yes, could do it based on individual likelihood and numeric gradient)
summary(b)
