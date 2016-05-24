library(kappalab)

x <- runif(14)
for (i in 2:14)
    x[i] <- x[i] + x[i-1]
mu <- capacity(c(0,0,x))
a <- Mobius(mu)
stopifnot(abs(entropy(mu) - entropy(a)) < 1e-6)
mu <- lower.capacity(4)
stopifnot(abs(entropy(mu) - entropy(as.capacity(as.set.func(mu)))) < 1e-6)
mu <- uniform.capacity(4)
stopifnot(abs(entropy(mu) - entropy(as.capacity(as.set.func(mu)))) < 1e-6)




