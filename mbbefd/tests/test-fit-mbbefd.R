library(mbbefd)
library(fitdistrplus)

n <- 1e3
nboot <- 1000
nboot <- 10
set.seed(123456)
lossrate <- rmbbefd(n, 1/2, 1/10)

f1 <- fitDR(lossrate, "mbbefd")

summary(f1)
cdfcomp(f1, do.points=FALSE)
qqcomp(f1)


# llsurface(plot.min=c(0, 0), plot.max=c(2, 1/2), plot.arg=c("a", "b"), obs=lossrate, distr="mbbefd", nlevels=25)
# points(f1$estimate["a"], f1$estimate["b"], pch="+", col="red")
# points(1/2, 1/10, pch="x", col="black")

b1 <- bootDR(f1, niter=nboot, silent=TRUE)
plot(b1, enhance=TRUE, trueval=c(1/2, 1/10))


f2 <- fitDR(lossrate, "mbbefd", method="tlmme")
summary(f2)



set.seed(123456)
lossrate <- rmbbefd(n, -1/2, 5)


f1 <- fitDR(lossrate, "mbbefd")
summary(f1)

b1 <- bootDR(f1, niter=nboot, silent=TRUE)
plot(b1, enhance=TRUE, trueval=c(-1/2, 5))

