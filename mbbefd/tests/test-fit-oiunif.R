library(mbbefd)
library(fitdistrplus)


#oiunif
n <- 1e3
nboot <- 1000
nboot <- 10
x <- roiunif(n, 1/6)
f1 <- fitDR(x, "oiunif", method="mle")
summary(f1)
summary(fitdist(x, "oiunif", method="mle", start=list(p1=1/2))) #check

b1 <- bootDR(f1, niter=nboot, silent=FALSE)
summary(b1)

plot(b1)
abline(v=1/6, col="red")

hist(b1$estim[,1])
abline(v=1/6, col="red")


f2 <- fitDR(x, "oiunif", method="tlmme")


