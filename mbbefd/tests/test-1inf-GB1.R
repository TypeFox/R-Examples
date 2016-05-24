library(mbbefd)

#test of GB1 distribution

#integral of the improper density
integrate(doigbeta, 0, 1, shape0=1, shape1=3, shape2=3/2, p1=1/3)
integrate(doigbeta, 0, 1, shape0=1/2, shape1=3, shape2=3/2, p1=1/3)
integrate(doigbeta, 0, 1, shape0=2, shape1=3, shape2=3/2, p1=2/3)


#RNG
n <- 1e4
x <- roigbeta(n, shape0=2, shape1=3, shape2=3/2, p1=1/3)
y <- roigbeta(n, shape0=pi, shape1=3, shape2=3/2, p1=2/3)

c(etl(x), tloigbeta(shape0=2, shape1=3, shape2=3/2, p1=1/3))
c(etl(y), tloigbeta(shape0=pi, shape1=3, shape2=3/2, p1=2/3))

#test CDF
z <- 0:10/10
cbind(ecdf(x)(z), poigbeta(z, shape0=2, shape1=3, shape2=3/2, p1=1/3))
cbind(ecdf(y)(z), poigbeta(z, shape0=pi, shape1=3, shape2=3/2, p1=2/3))


#mean
c(mean(x), moigbeta(1, shape0=2, shape1=3, shape2=3/2, p1=1/3))
c(mean(y), moigbeta(1, shape0=pi, shape1=3, shape2=3/2, p1=2/3))

#raw moment
for(i in 2:4)
{
  cat("E(X^", i, ")\n", sep="")
  print(c(mean(x^i), moigbeta(i, shape0=2, shape1=3, shape2=3/2, p1=1/3)))
  print(c(mean(y^i), moigbeta(i, shape0=pi, shape1=3, shape2=3/2, p1=2/3)))
}


#test EC
cbind(eecf(x)(z), ecoigbeta(z, shape0=2, shape1=3, shape2=3/2, p1=1/3))

cbind(eecf(y)(z), ecoigbeta(z, shape0=pi, shape1=3, shape2=3/2, p1=2/3))




