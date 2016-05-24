library(sensR)

## Testing Tetrad link function:
d.primePwr(1, sample.size = 50, method = "tetrad")

d.primeSS(1, target.power = 0.90, method = "tetrad")
discrim(10, 15, method = "tetrad", statistic = "score")

set.seed(12345)
discrimSim(sample.size = 10, replicates = 3, d.prime = 2,
           method = "tetrad", sd.indiv = 1)

rescale(pd = c(.6,.7), std.err = c(.2, NA), method="tetrad")
psyfun(2, method = "tetrad")
psyinv(0.8, method = "tetrad")
psyderiv(2, method = "tetrad")

x <- c(3,2,6,8,3,4,6,0,9,9,0,2,1,2,8,9,5,7)
n <- c(10,9,8,9,8,6,9,10,10,10,9,9,10,10,10,10,9,10)
dat <- data.frame(x, n)

(bb <- betabin(dat, method = "tetrad"))
(bb <- betabin(dat, corrected = FALSE, method = "tetrad"))
summary(bb)
vcov(bb)
logLik(bb)
AIC(bb)
coef(bb)

## Testing that the A-not A confidence interval finds the right
## confint method for glm objects in MASS:
m1 <- AnotA(10, 20, 3, 20)

## Make ci.res:
## ci.res := dput(as.vector(confint(m1)))
## Compare with current results:
ci.res <- c(-0.550522667455557, 0.190017579324113, 0.550522667455557,
            1.93313337684111)
stopifnot(isTRUE(all.equal(ci.res, as.vector(confint(m1)))))
