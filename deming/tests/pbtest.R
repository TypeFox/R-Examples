library(deming)

# Data set from Sen.  Since it is monotone the Thiel-Sen and PB regressions
#  are identical
#
sen <- data.frame(t=c(1,  2,  3,  4, 10, 12, 18),
                  y=c(9, 15, 19, 20, 45, 55, 78))
fit <- pbreg(y ~ t, data=sen, conf=.93)
all.equal(coef(fit), c("(Intercept)"=6, "t"=4))
all.equal(fit$ci[2,], c(26/7, 35/8), check.attributes=FALSE)

fit2 <- pbreg(y ~ t, data=sen, conf=.93, method=2)
fit3 <- pbreg(y ~ t, data=sen, conf=.93, method=3)
all.equal(coef(fit2), coef(fit))
all.equal(coef(fit3), coef(fit))
all.equal(fit2$ci, fit$ci)
all.equal(fit3$ci, fit$ci)

# Simple data set with noise
pdata <- data.frame(x=1:6, y=c(2,1.2, 4.1,3.5,6.3, 3))
fit <- pbreg(y~x, pdata)
all.equal(coef(fit), c(1/8, 1.05), check.attributes=FALSE)

fit3 <- pbreg(y ~ x, pdata, method=3)
all.equal(coef(fit3), coef(fit))

fit2 <- pbreg(y ~ x, pdata, method=2)
temp <- tan(mean(c(atan(2.2/2), atan(4.3/4))))
all.equal(as.vector(coef(fit2)), c(median(pdata$y- pdata$x*temp), temp))
