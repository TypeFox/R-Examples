library(kappalab)

mu <- card.set.func(rnorm(10))
stopifnot(abs(mu@data - zeta(Mobius(mu))@data) < 1e-6)

mu <- card.set.func(runif(9))
stopifnot(abs(mu@data - as.card.set.func(as.set.func(mu))@data) < 1e-6)



