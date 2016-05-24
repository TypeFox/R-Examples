### R code from vignette source 'count-insolvency.Rnw'

###################################################
### code chunk number 1: count-insolvency.Rnw:11-14 (eval = FALSE)
###################################################
## library(catdata)
## data(insolvency)
## attach(insolvency)


###################################################
### code chunk number 2: count-insolvency.Rnw:19-22 (eval = FALSE)
###################################################
## ins1 <- glm(insolv ~ case + I(case^2), family=poisson(link=log), data=insolvency)
## summary(ins1)
## # plot(ins1)


###################################################
### code chunk number 3: count-insolvency.Rnw:25-27 (eval = FALSE)
###################################################
## plot(case, insolv)
## points(ins1$fitted.values, type="l")


###################################################
### code chunk number 4: count-insolvency.Rnw:31-34 (eval = FALSE)
###################################################
## ins2 <- glm(insolv ~ case + I(case^2), family=quasipoisson, data=insolvency)
## summary(ins2)
## # plot(ins2)


###################################################
### code chunk number 5: count-insolvency.Rnw:37-40 (eval = FALSE)
###################################################
## library(MASS)
## ins3 <- glm.nb(insolv ~ case + I(case^2),data=insolvency)
## summary(ins3)


###################################################
### code chunk number 6: count-insolvency.Rnw:43-45 (eval = FALSE)
###################################################
## ins4 <- glm(insolv ~ case + I(case^2), family=gaussian(link=log), data=insolvency)
## summary(ins4)


