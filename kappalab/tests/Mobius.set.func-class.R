library(kappalab)

a <- Mobius.set.func(rnorm(32),5,5)
mu <- zeta(a)
stopifnot(abs(mu@data - zeta(Mobius(mu))@data) < 1e-6)

stopifnot(abs(Shapley.value(a) - Shapley.value(mu)) < 1e-6)
ii1 <- interaction.indices(a)
ii2 <- interaction.indices(mu)
ii1[is.na(ii1)] <- 1
ii2[is.na(ii2)] <- 1
stopifnot(abs(ii1 -ii2) < 1e-6)

mu <- set.func(c(0,1,1,1,2,2,2,3))
a <- Mobius(mu)
stopifnot(abs(a@data - as.Mobius.set.func(as.Mobius.card.set.func(a))@data) < 1e-6)


