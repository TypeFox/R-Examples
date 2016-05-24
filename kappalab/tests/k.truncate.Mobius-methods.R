library(kappalab)

## a set function
mu <- set.func(c(0,1,1,1,2,2,2,3))

## 1-truncate it
a <- k.truncate.Mobius(mu,1)
stopifnot(abs(zeta(a)@data - mu@data) < 1e-6)
## 2-truncate it
a <- k.truncate.Mobius(Mobius(mu),1)
stopifnot(abs(zeta(a)@data - mu@data) < 1e-6)




