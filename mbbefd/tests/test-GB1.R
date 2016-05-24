library(mbbefd)

#test of GB1 distribution




#integral of the density
integrate(dgbeta, 0, 1, shape0=1, shape1=3, shape2=3/2)
integrate(dgbeta, 0, 1, shape0=1/2, shape1=3, shape2=3/2)
integrate(dgbeta, 0, 1, shape0=2, shape1=3, shape2=3/2)

z <- 0:10/10
cbind(pi*z^(3*pi-1)*(1-z^pi)^(3/2-1)/beta(3, 3/2), dgbeta(z, pi, 3, 3/2))
cbind(dbeta(z^pi, 3, 3/2)*pi*z^(pi-1), dgbeta(z, pi, 3, 3/2))
cbind(pbeta(z^pi, 3, 3/2), pgbeta(z, pi, 3, 3/2))
cbind(qbeta(z, 3, 3/2)^(1/pi), qgbeta(z, pi, 3, 3/2))


#log argument

cbind(dbeta(z, 3, 3/2, log=TRUE), dgbeta(z, 1, 3, 3/2, log=TRUE))
cbind(log(dgbeta(z, pi, 3, 3/2, log=FALSE)), dgbeta(z, pi, 3, 3/2, log=TRUE))
cbind(log(pgbeta(z, pi, 3, 3/2, log=FALSE)), pgbeta(z, pi, 3, 3/2, log=TRUE))
cbind(log(pgbeta(z, pi, 3, 3/2, log=FALSE, lower=TRUE)), pgbeta(z, pi, 3, 3/2, log=TRUE, lower=TRUE))
cbind(qgbeta(z, pi, 3, 3/2, log=FALSE), qgbeta(log(z), pi, 3, 3/2, log=TRUE))
cbind(qgbeta(z, pi, 3, 3/2, log=FALSE, lower=TRUE), qgbeta(log(z), pi, 3, 3/2, log=TRUE, lower=TRUE))

#RNG
n <- 1e4
x <- rgbeta(n, shape0=2, shape1=3, shape2=3/2)
y <- rgbeta(n, shape0=pi, shape1=3, shape2=3/2)

#test density
z <- 0:100/100

plot(density(x)); lines(z, dgbeta(z, shape0=2, shape1=3, shape2=3/2), col="red")
plot(density(y)); lines(z, dgbeta(z, shape0=pi, shape1=3, shape2=3/2), col="red")

#mode
modeGB1 <- function(shape0, shape1, shape2)
{
  if(shape1+shape2-1/shape0>1)
    ((shape1-1/shape0)/(shape1+shape2-1/shape0-1))^(1/shape0)
  else
    NaN
}
c(modeGB1(2, 3, 3/2), density(x)$x[which.max(density(x)$y)])
c(modeGB1(pi, 3, 3/2), density(y)$x[which.max(density(y)$y)])

#test CDF
z <- 0:10/10
cbind(ecdf(x)(z), pgbeta(z, shape0=2, shape1=3, shape2=3/2))

cbind(ecdf(y)(z), pgbeta(z, shape0=pi, shape1=3, shape2=3/2))

#plot(ecdf(x)); lines(z, pgbeta(z, shape0=2, shape1=3, shape2=3/2), col="red")
#plot(ecdf(y)); lines(z, pgbeta(z, shape0=pi, shape1=3, shape2=3/2), col="red")

#mean
c(mean(x), mgbeta(1, shape0=2, shape1=3, shape2=3/2))
c(mean(y), mgbeta(1, shape0=pi, shape1=3, shape2=3/2))

#raw moment
for(i in 2:4)
{
  cat("E(X^", i, ")\n", sep="")
  print(c(mean(x^i), mgbeta(i, shape0=2, shape1=3, shape2=3/2)))
  print(c(mean(y^i), mgbeta(i, shape0=pi, shape1=3, shape2=3/2)))
}

#test limited expected value
d <- 1/2
s0 <- 2
s1 <- 3
s2 <- 3/2

mean(pmin(x, d))

f <- function(x, d, shape0) dgbeta(x, shape0=shape0, shape1=3, shape2=3/2)*pmin(x,d)
integrate(f, 0, 1, d=d, shape0=s0)



#test EC
f <- function(x, d, shape0) dgbeta(x, shape0=shape0, shape1=3, shape2=3/2)*pmin(x, d)/mgbeta(1, shape0=shape0, shape1=3, shape2=3/2)
F <- function(d, shape0) sapply(d, function(d) integrate(f, 0, 1, d=d, shape0=shape0)$value)

cbind(eecf(x)(z), ecgbeta(z, shape0=2, shape1=3, shape2=3/2), F(z, shape0=2))

cbind(eecf(y)(z), ecgbeta(z, shape0=pi, shape1=3, shape2=3/2), F(z, shape0=pi))




