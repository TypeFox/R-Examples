
library(rockchalk)

## Replicate some R classics.  The budworm.lg data from predict.glm
## will work properly after re-formatting the information as a data.frame:

## example from Venables and Ripley (2002, pp. 190-2.)
df <- data.frame(ldose = rep(0:5, 2),
                 sex = factor(rep(c("M", "F"), c(6, 6))),
                 SF.numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16))
df$SF.numalive = 20 - df$SF.numdead

budworm.lg <- glm(cbind(SF.numdead, SF.numalive) ~ sex*ldose, data = df,
                  family = binomial)


plotCurves(budworm.lg, plotx = "ldose", modx = "sex", interval = "confidence",
           ylim = c(0, 1))

## See infert
model2 <- glm(case ~ age + parity + education + spontaneous + induced,
              data = infert, family = binomial())

## Curvature so slight we can barely see it
model2pc1 <- plotCurves(model2, plotx = "age", modx = "education",
                        interval = "confidence", ylim = c(0, 1))


model2pc2 <- plotCurves(model2, plotx = "age", modx = "education",
                        modxVals = levels(infert$education)[1],
                        interval = "confidence", ylim = c(0, 1))


model2pc2 <- plotCurves(model2, plotx = "age", modx = "education",
                        modxVals = levels(infert$education)[c(2,3)],
                        interval = "confidence", ylim = c(0, 1))

model2pc2 <- plotCurves(model2, plotx = "age", modx = "education",
                        modxVals = levels(infert$education)[c(2,3)],
                        ylim = c(0, 1), type = "response")





## Manufacture some data
set.seed(12345)
N <- 500
dat <- genCorrelatedData2(N = 500, means = c(5, 0, 0, 0), sds = rep(1, 4),
                          rho = 0.2, beta = rep(1, 5), stde = 20)

dat$xcat1 <- gl(2, N/2, labels = c("Monster", "Human"))
dat$xcat2 <- cut(rnorm(N), breaks = c(-Inf, 0, 0.4, 0.9, 1, Inf),
             labels = c("R", "M", "D", "P", "G"))

###The design matrix for categorical variables, xcat numeric
dat$xcat1n <- with(dat, contrasts(xcat1)[xcat1, ])
dat$xcat2n <- with(dat, contrasts(xcat2)[xcat2, ])


stde <- 2
dat$y <- with(dat, 0.03 + 11.5 * log(x1) * contrasts(dat$xcat1)[dat$xcat1] +
              0.1 * x2 + 0.04 * x2^2 + stde*rnorm(N))

stde <- 1
dat$y2 <- with(dat, 0.03 + 0.1 * x1 + 0.1 * x2 + 0.25 * x1 * x2 + 0.4 * x3 -
               0.1 * x4 + stde * rnorm(N))
stde <- 8
dat$y3 <- with(dat, 3 + 0.5 * x1 + 1.2 * (as.numeric(xcat1)-1) +
-0.8 * (as.numeric(xcat1)-1) * x1 +  stde * rnorm(N))

stde <- 8
dat$y4 <- with(dat, 3 + 0.5 * x1 +
               contrasts(dat$xcat2)[dat$xcat2, ] %*% c(0.1, -0.2, 0.3, 0.05)  +
               stde * rnorm(N))


## Curvature with interaction
m1 <- lm(y ~ log(x1)*xcat1 + x2 + I(x2^2), data=dat)
summary(m1)

## First, with no moderator
plotCurves(m1, plotx = "x1")

plotCurves(m1, plotx = "x1", modx = "xcat1")

## ## Verify that plot by comparing against a manually contructed alternative
## par(mfrow=c(1,2))
## plotCurves(m1, plotx = "x1", modx = "xcat1")
## newdat <- with(dat, expand.grid(x1 = plotSeq(x1, 30), xcat1 = levels(xcat1)))
## newdat$x2 <-  with(dat, mean(x2, na.rm = TRUE))
## newdat$m1p <- predict(m1, newdata = newdat)
## plot( y ~ x1, data = dat, type = "n", ylim = magRange(dat$y, c(1, 1.2)))
## points( y ~ x1, data = dat, col = dat$xcat1, cex = 0.4, lwd = 0.5)
## by(newdat, newdat$xcat1, function(dd) {lines(dd$x1, dd$m1p)})
## legend("topleft", legend=levels(dat$xcat1), col = as.numeric(dat$xcat1), lty = 1)
## par(mfrow = c(1,1))
## ##Close enough!


plotCurves(m1, plotx = "x2", modx = "x1")

plotCurves(m1, plotx = "x2", modx = "xcat1")

plotCurves(m1, plotx = "x2", modx = "xcat1", interval = "conf")


m2 <- lm(y ~ log(x1)*xcat1 + xcat1*(x2 + I(x2^2)), data = dat)
summary(m2)
plotCurves(m2, plotx = "x2", modx = "xcat1")

plotCurves(m2, plotx  ="x2", modx = "x1")


m3a <- lm(y ~ poly(x2, 2) + xcat1, data = dat)

plotCurves(m3a, plotx = "x2")
plotCurves(m3a, plotx = "x2", modx = "xcat1")
#OK

m4 <- lm(log(y+10) ~ poly(x2, 2)*xcat1 + x1, data = dat)
summary(m4)
plotCurves(m4, plotx = "x2")

plotCurves(m4, plotx  ="x2", modx = "xcat1")

plotCurves(m4, plotx = "x2", modx = "x1")

plotCurves(m4, plotx = "x2", modx = "xcat1")

plotCurves(m4, plotx = "x2", modx = "xcat1", modxVals = c("Monster"))


##ordinary interaction
m5 <- lm(y2 ~ x1*x2 + x3 +x4, data = dat)
summary(m5)
plotCurves(m5, plotx = "x1", modx = "x2")
plotCurves(m5, plotx = "x1", modx = "x2", modxVals = c( -2, -1, 0, 1, 2))
plotCurves(m5, plotx = "x1", modx = "x2", modxVals = c(-2))
plotCurves(m5, plotx = "x1", modx = "x2", modxVals = "std.dev.")
plotCurves(m5, plotx = "x1", modx = "x2", modxVals = "quantile")
plotCurves(m5, plotx = "x3", modx = "x2")


library(car)
mc1 <- lm(statusquo ~ income * sex, data = Chile)
summary(mc1)
plotCurves(mc1, plotx = "income")
plotCurves(mc1, modx = "sex", plotx = "income")
plotCurves(mc1, modx = "sex", plotx = "income", modxVals = "M")

mc2 <- lm(statusquo ~ region * income, data =  Chile)
summary(mc2)
plotCurves(mc2, modx = "region", plotx = "income")
plotCurves(mc2, modx = "region", plotx = "income",
           modxVals = levels(Chile$region)[c(1,4)])
plotCurves(mc2, modx = "region", plotx = "income", modxVals = c("S","M","SA"))
plotCurves(mc2, modx = "region", plotx = "income", modxVals = c("S","M","SA"),
           interval = "conf")

plotCurves(mc2, modx = "region", plotx = "income", plotPoints = FALSE)


mc3 <- lm(statusquo ~ region * income + sex + age, data =  Chile)
summary(mc3)
plotCurves(mc3, modx = "region", plotx = "income")


mc4 <- lm(statusquo ~ income * (age + I(age^2)) + education + sex + age, data = Chile)
summary(mc4)
plotCurves(mc4, plotx = "age")
plotCurves(mc4, plotx = "age", interval = "conf")


plotCurves(mc4, plotx = "age", modx = "income")
plotCurves(mc4, plotx = "age", modx = "income", plotPoints = FALSE)

plotCurves(mc4,  plotx = "income", modx = "age")
plotCurves(mc4,  plotx = "income", modx = "age", n = 8)

plotCurves(mc4,  plotx = "income", modx = "age", modxVals = "std.dev.")
plotCurves(mc4, modx = "income", plotx = "age", plotPoints = FALSE)
