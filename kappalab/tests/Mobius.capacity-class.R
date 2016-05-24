library(kappalab)

x <- runif(15)
for (i in 2:15)
    x[i] <- x[i] + x[i-1]
mu <- capacity(c(0,x))
a <- Mobius(mu)

stopifnot(abs(orness(a) - orness(mu)) < 1e-6)
stopifnot(abs(favor(a) - favor(mu)) < 1e-6)
stopifnot(abs(veto(a) - veto(mu)) < 1e-6)
stopifnot(abs(variance(a) - variance(mu)) < 1e-6)
stopifnot(abs(variance(a) - variance(mu)) < 1e-6)
stopifnot(abs(mean(favor(mu)) - orness(mu)) < 1e-6)
stopifnot(abs(mean(veto(mu)) - 1 + orness(mu)) < 1e-6)
stopifnot(abs(mean(favor(a)) - orness(a)) < 1e-6)
stopifnot(abs(mean(veto(a)) - 1 + orness(a)) < 1e-6)



