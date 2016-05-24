library(kappalab)

x <- runif(15)
for (i in 2:15)
    x[i] <- x[i] + x[i-1]
mu <- capacity(c(0,x))
a <- Mobius(mu)
stopifnot(abs(variance(mu) - variance(a)) < 1e-6)

mu <- card.capacity(c(0,x[1:6]))
stopifnot(abs(variance(mu) - variance(as.capacity(as.set.func(mu)))) < 1e-6)

