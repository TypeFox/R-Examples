library(kappalab)

x <- runif(15)
for (i in 2:15)
    x[i] <- x[i] + x[i-1]
mu <- capacity(c(0,x))
a <- Mobius(mu)
stopifnot(abs(mu@data - zeta(Mobius(mu))@data) < 1e-6)

mu <- game(c(0,runif(31)))
stopifnot(abs(mu@data - zeta(Mobius(mu))@data) < 1e-6)

mu <- set.func(rnorm(64))
stopifnot(abs(mu@data - zeta(Mobius(mu))@data) < 1e-6)

mu <- card.set.func(rnorm(7,10,2))
stopifnot(abs(mu@data - zeta(Mobius(mu))@data) < 1e-6)
stopifnot(abs(mu@data - as.card.set.func(as.set.func(mu))@data) < 1e-6)


