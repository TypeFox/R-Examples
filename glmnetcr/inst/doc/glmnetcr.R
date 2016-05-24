### R code from vignette source 'glmnetcr.Rnw'

###################################################
### code chunk number 1: glmnetcr.Rnw:74-75
###################################################
options(prompt="R> ", continue="+  ", width=70, useFancyQuotes=FALSE)


###################################################
### code chunk number 2: glmnetcr.Rnw:77-85
###################################################
library("glmnetcr")
data("diabetes")
dim(diabetes)
names(diabetes)[1:10]
summary(diabetes$y)
x <- diabetes[, 2:dim(diabetes)[2]]
y <- diabetes$y
fit <- glmnet.cr(x, y)


###################################################
### code chunk number 3: glmnetcr.Rnw:88-89
###################################################
print(fit)


###################################################
### code chunk number 4: glmnetcr.Rnw:93-94
###################################################
plot(fit, xvar = "step", type = "bic")


###################################################
### code chunk number 5: glmnetcr.Rnw:102-103
###################################################
plot(fit, xvar = "step", type = "coefficients")


###################################################
### code chunk number 6: glmnetcr.Rnw:110-112
###################################################
plot(fit, xvar = "step", type = "bic")
plot(fit, xvar = "step", type = "coefficients")


###################################################
### code chunk number 7: glmnetcr.Rnw:116-120
###################################################
BIC.step <- select.glmnet.cr(fit)
BIC.step
AIC.step <- select.glmnet.cr(fit, which = "AIC")
AIC.step


###################################################
### code chunk number 8: glmnetcr.Rnw:123-127
###################################################
coefficients<-coef(fit, s = BIC.step)
coefficients$a0
sum(coefficients$beta != 0)
nonzero.glmnet.cr(fit, s = BIC.step)


###################################################
### code chunk number 9: glmnetcr.Rnw:131-132
###################################################
fit <- glmnet.cr(x, y, method = "forward")


###################################################
### code chunk number 10: glmnetcr.Rnw:135-137
###################################################
BIC.step <- select.glmnet.cr(fit)
BIC.step


###################################################
### code chunk number 11: glmnetcr.Rnw:141-145
###################################################
coefficients<-coef(fit, s = BIC.step)
coefficients$a0
sum(coefficients$beta != 0)
nonzero.glmnet.cr(fit, s = BIC.step)


###################################################
### code chunk number 12: glmnetcr.Rnw:150-153
###################################################
hat<-fitted(fit, s = BIC.step)
names(hat)
table(hat$class, y)


