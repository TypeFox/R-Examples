library(kappalab)

x <- runif(15)
for (i in 2:15)
    x[i] <- x[i] + x[i-1]
mu <- capacity(c(0,x))
a <- Mobius(mu)
f <- c(0.1,0.9,0.3,0.8)
stopifnot(abs(Choquet.integral(mu,f) - Choquet.integral(a,f)) < 1e-6)



