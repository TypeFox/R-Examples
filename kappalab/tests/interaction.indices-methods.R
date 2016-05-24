library(kappalab)

mu <- set.func(rnorm(32))
a <- Mobius(mu)
ii1 <- interaction.indices(a)
ii2 <- interaction.indices(mu)
ii1[is.na(ii1)] <- 1
ii2[is.na(ii2)] <- 1
stopifnot(abs(ii1 -ii2) < 1e-6)

mu <- upper.capacity(6)
stopifnot(abs(interaction.indices(mu)[1,1] - interaction.indices(as.set.func(mu))[1,2]) < 1e-6)

mu <- card.set.func(rnorm(12))
mu2 <- as.set.func(mu)
ii <- interaction.indices(mu2)
stopifnot(abs(ii[1,2] -interaction.indices(mu)) < 1e-6)
