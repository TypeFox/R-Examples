library(mbbefd)

#test of MBBEFD(a,b) distribution
n <- 1e4

set.seed(567)
x <- rmbbefd(n, 2, 1/2)
y <- rmbbefd(n, -1/2, 2)



#test CDF
z <- 0:8/8
cbind(ecdf(x)(z), pmbbefd(z, 2, 1/2))

cbind(ecdf(y)(z), pmbbefd(z, -1/2, 2))

#test EC
cbind(eecf(x)(z), ecmbbefd(z, 2, 1/2))

cbind(eecf(y)(z), ecmbbefd(z, -1/2, 2))

#test mean
mean(x)
mmbbefd(1, 2, 1/2)

mean(y)
mmbbefd(1, -1/2, 2)

#total loss
etl(x)
tlmbbefd(2, 1/2)


etl(y)
tlmbbefd(-1/2, 2)


#test quantile
cbind(quantile(y, probs=0:10/10), qmbbefd(0:10/10, -1/2, 2))
qmbbefd(1/2, -1/2, 2)


z <- seq(0, 1, length=101)
plot(z, pmbbefd(z, -1/2, 2), type="l", ylim=c(0, 1-tlmbbefd(-1/2, 2)))

plot(z, qmbbefd(z, -1/2, 2), type="l", xlim=c(0, 1-tlmbbefd(-1/2, 2)))



#test density

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(x), ylim=c(0,1))
lines(z, dmbbefd(z, 2, 1/2), col="red")


plot(density(y), ylim=c(0,1))
lines(z, dmbbefd(z, -1/2, 2), col="red")
