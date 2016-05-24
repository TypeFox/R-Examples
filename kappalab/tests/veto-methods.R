library(kappalab)

x <- runif(15)
for (i in 2:15)
    x[i] <- x[i] + x[i-1]
mu <- capacity(c(0,x))
a <- Mobius(mu)
stopifnot(abs(favor(mu) - favor(a)) < 1e-6)

mu <- card.capacity(c(0,x[1:7]))
stopifnot(abs(favor(mu) - favor(as.capacity(mu))) < 1e-6)



