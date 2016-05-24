library(mbbefd)

#test of shifted truncated pareto distribution
n <- 1e4

x <- roiunif(n, p1=1/2)
y <- roiunif(n, p1=1/3)

#test CDF
z <- 0:10/10
cbind(ecdf(x)(z), poiunif(z, 1/2))

cbind(ecdf(y)(z), poiunif(z, 1/3))

#total loss
c(etl(x), tloiunif(1/2))
c(etl(y), tloiunif(1/3))

#mean
c(mean(x), moiunif(1, 1/2))
c(mean(y), moiunif(1, 1/3))


#test EC
cbind(eecf(x)(z), ecoiunif(z, 1/2))
cbind(eecf(y)(z), ecoiunif(z, 1/3))



#plots
n <- 1e2
x <- roiunif(n, p1=1/2)
y <- roiunif(n, p1=1/3)


plot(eecf(x), do.points=FALSE)
v <- seq(0, 1, length=101)
lines(v, ecoiunif(v, 1/2), lty=3, col="red")


plot(eecf(y), do.points=FALSE)
v <- seq(0, 1, length=101)
lines(v, ecoiunif(v, 1/3), lty=3, col="red")

