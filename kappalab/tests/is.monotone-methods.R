library(kappalab)

## a monotone set function
mu <- set.func(c(0,1,1,1,2,2,2,3))
mu
is.monotone(mu)

## the Mobius representation of a monotone set function
a <- Mobius.set.func(c(0,1,2,1,3,1,2,1,2,3,1),4,2)
is.monotone(a)

## non-monotone examples
mu <- set.func(c(0,-7:7))
is.monotone(mu,verbose=TRUE)
a <- Mobius(mu)
is.monotone(a,verbose=TRUE)



