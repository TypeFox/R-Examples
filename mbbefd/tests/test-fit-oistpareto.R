library(mbbefd)
library(fitdistrplus)


#oistpareto
n <- 1e3
nboot <- 1000
nboot <- 10
set.seed(12345)
x <- roistpareto(n, 2, 1/6)

f1 <- fitDR(x, "oistpareto", method="mle")
summary(f1)
summary(fitdist(x, "oistpareto", method="mle", start=list(a=1/mean(x), p1=etl(x))))#check

b1 <- bootDR(f1, niter=nboot)
summary(b1)

plot(b1, enhance=TRUE, trueval=c(2, 1/6))

hist(b1$estim[,1])
hist(b1$estim[,2])

f2 <- fitDR(x, "oistpareto", method="tlmme")
summary(f2)

gofstat(list(f1, f2))
cdfcomp(list(f1, f2), do.points=FALSE)
ppcomp(list(f1, f2))
