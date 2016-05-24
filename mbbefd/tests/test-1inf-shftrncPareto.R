library(mbbefd)

#test of shifted truncated pareto distribution
n <- 1e4

x <- roistpareto(n, 2, p1=1/2)
y <- roistpareto(n, 1/2, p1=1/3)

#test CDF
z <- 0:10/10
cbind(ecdf(x)(z), poistpareto(z, 2, 1/2))

cbind(ecdf(y)(z), poistpareto(z, 1/2, 1/3))

#total loss
c(etl(x), tloistpareto(2, 1/2))
c(etl(y), tloistpareto(1/2, 1/3))

#mean
c(mean(x), moistpareto(1, 2, 1/2))
c(mean(y), moistpareto(1, 1/2, 1/3))


#test EC
cbind(eecf(x)(z), ecoistpareto(z, 2, 1/2))
cbind(eecf(y)(z), ecoistpareto(z, 1/2, 1/3))



#plots
n <- 1e2
x <- roistpareto(n, 2, p1=1/2)
y <- roistpareto(n, 1/2, p1=1/3)


plot(eecf(x), do.points=FALSE)
v <- seq(0, 1, length=101)
lines(v, ecoistpareto(v, 2, 1/2), lty=3, col="red")


plot(eecf(y), do.points=FALSE)
v <- seq(0, 1, length=101)
lines(v, ecoistpareto(v, 1/2, 1/3), lty=3, col="red")

