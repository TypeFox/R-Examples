### R code from vignette source 'MOB.Rnw'

###################################################
### code chunk number 1: setup
###################################################
require("party")
options(useFancyQuotes = FALSE)


###################################################
### code chunk number 2: MOB.Rnw:232-233
###################################################
library("party")


###################################################
### code chunk number 3: MOB.Rnw:238-239
###################################################
data("BostonHousing", package = "mlbench")


###################################################
### code chunk number 4: MOB.Rnw:253-255
###################################################
BostonHousing$lstat <- log(BostonHousing$lstat)
BostonHousing$rm <- BostonHousing$rm^2


###################################################
### code chunk number 5: MOB.Rnw:270-272
###################################################
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, labels = c("no", "yes"))
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)


###################################################
### code chunk number 6: MOB.Rnw:283-285
###################################################
ctrl <- mob_control(alpha = 0.05, bonferroni = TRUE, minsplit = 40,
  objfun = deviance, verbose = TRUE)


###################################################
### code chunk number 7: MOB.Rnw:298-300
###################################################
fmBH <- mob(medv ~ lstat + rm | zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
  data = BostonHousing, control = ctrl, model = linearModel)


###################################################
### code chunk number 8: MOB.Rnw:309-310
###################################################
fmBH


###################################################
### code chunk number 9: MOB.Rnw:316-317 (eval = FALSE)
###################################################
## plot(fmBH)


###################################################
### code chunk number 10: BostonHousing-plot
###################################################
plot(fmBH)


###################################################
### code chunk number 11: MOB.Rnw:350-351
###################################################
coef(fmBH)


###################################################
### code chunk number 12: MOB.Rnw:357-358
###################################################
summary(fmBH, node = 7)


###################################################
### code chunk number 13: MOB.Rnw:365-366
###################################################
sctest(fmBH, node = 7)


###################################################
### code chunk number 14: MOB.Rnw:372-375
###################################################
mean(residuals(fmBH)^2)
logLik(fmBH)
AIC(fmBH)


###################################################
### code chunk number 15: MOB.Rnw:378-380
###################################################
nt <- NROW(coef(fmBH))
nk <- NCOL(coef(fmBH))


###################################################
### code chunk number 16: MOB.Rnw:396-398
###################################################
data("PimaIndiansDiabetes2", package = "mlbench")
PimaIndiansDiabetes <- na.omit(PimaIndiansDiabetes2[,-c(4, 5)])


###################################################
### code chunk number 17: MOB.Rnw:418-420
###################################################
fmPID <- mob(diabetes ~ glucose | pregnant + pressure + mass + pedigree + age,
  data = PimaIndiansDiabetes, model = glinearModel, family = binomial())


###################################################
### code chunk number 18: MOB.Rnw:425-426 (eval = FALSE)
###################################################
## plot(fmPID)


###################################################
### code chunk number 19: PimaIndiansDiabetes-plot
###################################################
plot(fmPID)


###################################################
### code chunk number 20: MOB.Rnw:460-462
###################################################
coef(fmPID)
exp(coef(fmPID)[,2])


###################################################
### code chunk number 21: MOB.Rnw:465-466
###################################################
risk <- round(100 * (exp(coef(fmPID)[,2])-1), digits = 1)


