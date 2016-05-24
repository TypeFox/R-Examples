library(kappalab)

x <- runif(15)
for (i in 2:15)
    x[i] <- x[i] + x[i-1]
mu <- set.func(c(0,x)/x[15])
stopifnot(abs(sum(Shapley.value(mu)) - 1) < 1e-6)

a <- Mobius(mu)
stopifnot(abs(Shapley.value(a)- Shapley.value(mu)) < 1e-6)
stopifnot(abs(sum(Shapley.value(a)) - 1) < 1e-6)

mu <- card.set.func(rnorm(11))
stopifnot(abs(Shapley.value(mu) - Shapley.value(as.set.func(mu))) < 1e-6)




