library(kappalab)

a <- Mobius.card.set.func(rnorm(13,3,4))
stopifnot(abs(a@data - Mobius(zeta(a))@data) < 1e-6)
stopifnot(abs(a@data - as.Mobius.card.set.func(as.Mobius.set.func(a))@data) < 1e-6)



