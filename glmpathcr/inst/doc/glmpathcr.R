### R code from vignette source 'glmpathcr.Rnw'

###################################################
### code chunk number 1: glmpathcr.Rnw:74-82
###################################################
library(glmpathcr)
data(diabetes)
dim(diabetes)
names(diabetes)[1:10]
summary(diabetes$y)
x <- diabetes[, 2:dim(diabetes)[2]]
y <- diabetes$y
fit <- glmpath.cr(x,y)


###################################################
### code chunk number 2: glmpathcr.Rnw:85-87
###################################################
summary(fit)
plot(fit, xvar = "step", type = "bic")


###################################################
### code chunk number 3: glmpathcr.Rnw:91-92
###################################################
plot(fit, xvar = "step", type = "bic")


###################################################
### code chunk number 4: glmpathcr.Rnw:98-102
###################################################
BIC.step <- model.select(fit)
BIC.step
AIC.step <- model.select(fit)
AIC.step


###################################################
### code chunk number 5: glmpathcr.Rnw:107-110
###################################################
coefficients<-coef(fit, s=BIC.step)
sum(coefficients!=0)
nonzero.coef(fit, s=BIC.step)


###################################################
### code chunk number 6: glmpathcr.Rnw:115-119
###################################################
pred <- predict(fit)
table(pred, y)
pred <- predict(fit, type="probs")
pred


###################################################
### code chunk number 7: glmpathcr.Rnw:124-125
###################################################
fit <- glmpath.cr(x, y, method="forward")


###################################################
### code chunk number 8: glmpathcr.Rnw:128-130
###################################################
coefficients<-coef(fit, s=BIC.step)
nonzero.coef(fit, s=BIC.step)


###################################################
### code chunk number 9: glmpathcr.Rnw:134-136
###################################################
pred <- predict(fit)
table(pred, y)


