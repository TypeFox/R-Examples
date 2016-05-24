library(kappalab)

mu <- capacity(c(0,rep(1,15)))
f <- c(0.1,0.1,0,0.9)
a <- Mobius(mu)

stopifnot(abs(Choquet.integral(mu,f)-Choquet.integral(a,f)) < 1e-6)
stopifnot(abs(orness(mu,f)-orness(a,f)) < 1e-6)

x <- runif(8)
for (i in 2:8)
    x[i] <- x[i] + x[i-1]
mu <- card.capacity(c(0,x))
f <- c(0.1,5.1,0,0.9,0.2,0.4,1.3,8)
stopifnot(abs(Choquet.integral(mu,f)-Choquet.integral(as.capacity(mu),f)) < 1e-6)
stopifnot(abs(orness(mu,f)-orness(as.capacity(mu),f)) < 1e-6)




