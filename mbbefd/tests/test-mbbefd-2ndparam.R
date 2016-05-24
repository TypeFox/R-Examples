library(mbbefd)

#test of MBBEFD(a,b) distribution
n <- 1e4

set.seed(567)
x <- rMBBEFD(n, 2, 1/5)
y <- rMBBEFD(n, 3, 2)



#test CDF
z <- 0:8/8
cbind(ecdf(x)(z), round(pMBBEFD(z, 2, 1/5), 6))

cbind(ecdf(y)(z), round(pMBBEFD(z, 3, 2), 6))

#test EC
cbind(eecf(x)(z), round(ecMBBEFD(z, 2, 1/5), 6))

cbind(eecf(y)(z), round(ecMBBEFD(z, 3, 2), 6))

#test mean
mean(x)
mMBBEFD(1, 2, 1/2)

mean(y)
mMBBEFD(1, 3, 2)

#total loss
etl(x)
tlMBBEFD(2, 1/2)


etl(y)
tlMBBEFD(3, 2)


#test quantile
cbind(quantile(x, probs=0:10/10), round(qMBBEFD(0:10/10, 2, 1/5), 6))
cbind(quantile(y, probs=0:10/10), round(qMBBEFD(0:10/10, 3, 2), 6))


z <- seq(0, 1, length=101)
plot(z, pMBBEFD(z, 3, 2), type="l", ylim=c(0, 1-tlMBBEFD(3, 2)))

plot(z, qMBBEFD(z, 3, 2), type="l", xlim=c(0, 1-tlMBBEFD(3, 2)))



#test density

integrate(dMBBEFD, 0, 1, g=2, b=1/5)
pMBBEFD(1-1e-6, 2, 1/5)


integrate(dMBBEFD, 0, 1, g=3, b=2)
pMBBEFD(1-1e-6, 3, 2)


z <- 0:8/8

d <- function(z) approxfun(density(x)$x, density(x)$y)(z)
cbind(d(z), dMBBEFD(z, 2, 1/5))

d <- function(z) approxfun(density(y)$x, density(y)$y)(z)
cbind(d(z), dMBBEFD(z, 3, 2))

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(x), ylim=c(0,1))
lines(z, dMBBEFD(z, 2, 1/5), col="red")


plot(density(y), ylim=c(0,3))
lines(z, dMBBEFD(z, 3, 2), col="red")
