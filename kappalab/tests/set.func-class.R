library(kappalab)

mu <- set.func(c(1:8,8:1)/8)
stopifnot(abs(mu@data - zeta(Mobius(mu))@data) < 1e-6)

a <- Mobius(mu)
stopifnot(abs(Shapley.value(a)- Shapley.value(mu)) < 1e-6)

mu <- set.func(c(0,1,1,1,2,2,2,3))
stopifnot(abs(mu@data - as.set.func(as.card.set.func(mu))@data) < 1e-6)
stopifnot(abs(Shapley.value(mu) - Shapley.value(as.card.set.func(mu))) < 1e-6)



