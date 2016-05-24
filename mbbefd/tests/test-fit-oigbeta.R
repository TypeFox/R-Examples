library(mbbefd)
library(fitdistrplus)


#oigbeta
n <- 1e3
nboot <- 1000
nboot <- 10
set.seed(12345)
x <- roigbeta(n, 3, 2, 5/2, 1/6)


f1 <- fitDR(x, "oigbeta", method="mle", control=list(trace=1, REPORT=1, maxit=500)) #
summary(f1)

b1 <- bootDR(f1, niter=nboot, silent=TRUE)
summary(b1)

plot(b1, enhance=TRUE, trueval=c(3, 2, 5/2, 1/6))

f2 <- fitDR(x, "oigbeta", method="tlmme")
summary(f2)

gofstat(list(f1, f2))
cdfcomp(list(f1, f2), do.points=FALSE, ylogscale = TRUE)
ppcomp(list(f1, f2))

