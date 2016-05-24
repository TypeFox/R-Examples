### R code from vignette source 'chemCal.Rnw'

###################################################
### code chunk number 1: chemCal.Rnw:38-42
###################################################
library(chemCal)
data(massart97ex3)
m0 <- lm(y ~ x, data = massart97ex3)
calplot(m0)


###################################################
### code chunk number 2: chemCal.Rnw:49-50
###################################################
plot(m0,which=3)


###################################################
### code chunk number 3: chemCal.Rnw:56-63
###################################################
attach(massart97ex3)
yx <- split(y, x)
ybar <- sapply(yx, mean)
s <- round(sapply(yx, sd), digits = 2)
w <- round(1 / (s^2), digits = 3)
weights <- w[factor(x)]
m <- lm(y ~ x, w = weights)


###################################################
### code chunk number 4: chemCal.Rnw:69-71
###################################################
inverse.predict(m, 15, ws=1.67)
inverse.predict(m, 90, ws = 0.145)


