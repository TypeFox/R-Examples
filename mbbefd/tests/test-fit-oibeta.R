library(mbbefd)
library(fitdistrplus)


#oibeta
n <- 1e3
nboot <- 1000
nboot <- 10
set.seed(12345)
x <- roibeta(n, 3, 2, 1/6)

f1 <- fitDR(x, "oibeta", method="mle")
summary(f1)

b1 <- bootDR(f1, niter=nboot)
summary(b1)

plot(b1, enhance=TRUE, trueval=c(3, 2, 1/6))

hist(b1$estim[,1])
hist(b1$estim[,2])
hist(b1$estim[,3])


f2 <- fitDR(x, "oigbeta", method="tlmme")
summary(f2)

gofstat(list(f1, f2))
cdfcomp(list(f1, f2), do.points=FALSE)
ppcomp(list(f1, f2))
