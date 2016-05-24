library(mbbefd)

#test of shifted truncated pareto distribution
n <- 1e4

x <- rstpareto(n, 2)
y <- rstpareto(n, 1/2)

#test CDF
z <- 0:10/10
cbind(ecdf(x)(z), pstpareto(z, 2))

cbind(ecdf(y)(z), pstpareto(z, 1/2))

#mean
c(mean(x), mstpareto(1, 2))
c(mean(y), mstpareto(1, 1/2))

#test EC
cbind(eecf(x)(z), ecstpareto(z, 2))

cbind(eecf(y)(z), ecstpareto(z, 1/2))


plot(eecf(x))
v <- seq(0, 1, length=101)
lines(v, ecstpareto(v, 2), lty=3, col="red")


plot(eecf(y))
v <- seq(0, 1, length=101)
lines(v, ecstpareto(v, 1/2), lty=3, col="red")

