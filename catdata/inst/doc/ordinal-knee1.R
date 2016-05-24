### R code from vignette source 'ordinal-knee1.Rnw'

###################################################
### code chunk number 1: ordinal-knee1.Rnw:12-13 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: ordinal-knee1.Rnw:19-23 (eval = FALSE)
###################################################
## rm(list=ls(all=TRUE))
## library(catdata)
## data(knee)
## attach(knee)


###################################################
### code chunk number 3: ordinal-knee1.Rnw:27-28 (eval = FALSE)
###################################################
## suppressWarnings(chisq.test(knee$Th,knee$R4))


###################################################
### code chunk number 4: ordinal-knee1.Rnw:32-34 (eval = FALSE)
###################################################
## Age <- Age - 30
## Age2 <- Age^2


###################################################
### code chunk number 5: ordinal-knee1.Rnw:38-41 (eval = FALSE)
###################################################
## R4 <- as.ordered(R4)
## Th <- as.factor(Th)
## Sex <- as.factor(Sex)


###################################################
### code chunk number 6: ordinal-knee1.Rnw:46-47 (eval = FALSE)
###################################################
## library(MASS)


###################################################
### code chunk number 7: ordinal-knee1.Rnw:52-53 (eval = FALSE)
###################################################
## polr1 <- polr(R4 ~ Th, method="logistic")


###################################################
### code chunk number 8: ordinal-knee1.Rnw:56-57 (eval = FALSE)
###################################################
## summary(polr1)


###################################################
### code chunk number 9: ordinal-knee1.Rnw:61-62 (eval = FALSE)
###################################################
## exp(-coef(polr1))


###################################################
### code chunk number 10: ordinal-knee1.Rnw:68-69 (eval = FALSE)
###################################################
## polr2 <- polr(R4 ~ Th + Sex + Age, method="logistic")


###################################################
### code chunk number 11: ordinal-knee1.Rnw:72-73 (eval = FALSE)
###################################################
## summary(polr2)


###################################################
### code chunk number 12: ordinal-knee1.Rnw:77-78 (eval = FALSE)
###################################################
## exp(-coef(polr2))


###################################################
### code chunk number 13: ordinal-knee1.Rnw:84-86 (eval = FALSE)
###################################################
## se <- summary(polr2)[1][[1]][1:3,2]
## (wald2 <- -coef(polr2)/se)


###################################################
### code chunk number 14: ordinal-knee1.Rnw:90-91 (eval = FALSE)
###################################################
## 1-pchisq(wald2^2, df=1)


###################################################
### code chunk number 15: ordinal-knee1.Rnw:96-97 (eval = FALSE)
###################################################
## polr3 <- update(polr2, ~. + Age2)


###################################################
### code chunk number 16: ordinal-knee1.Rnw:100-101 (eval = FALSE)
###################################################
## summary(polr3)


###################################################
### code chunk number 17: ordinal-knee1.Rnw:105-106 (eval = FALSE)
###################################################
## exp(-coef(polr3))


###################################################
### code chunk number 18: ordinal-knee1.Rnw:111-113 (eval = FALSE)
###################################################
## se <- summary(polr3)[1][[1]][1:4,2]
## (wald3 <- -coef(polr3)/se)


###################################################
### code chunk number 19: ordinal-knee1.Rnw:117-118 (eval = FALSE)
###################################################
## 1-pchisq(wald3^2, df=1)


###################################################
### code chunk number 20: ordinal-knee1.Rnw:129-130 (eval = FALSE)
###################################################
## library(VGAM)


###################################################
### code chunk number 21: ordinal-knee1.Rnw:133-136 (eval = FALSE)
###################################################
## m.vglm <- vglm(R4 ~ Th + Sex + Age + Age2, family = cumulative (link="logit", 
## parallel=TRUE))
## summary(m.vglm)


###################################################
### code chunk number 22: ordinal-knee1.Rnw:145-146 (eval = FALSE)
###################################################
## library(rms)


###################################################
### code chunk number 23: ordinal-knee1.Rnw:149-151 (eval = FALSE)
###################################################
## m.lrm <- lrm(R4 ~ Th + Sex + Age + Age2)
## m.lrm


###################################################
### code chunk number 24: ordinal-knee1.Rnw:156-157 (eval = FALSE)
###################################################
## detach(package:VGAM)


