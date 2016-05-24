### R code from vignette source 'qualint.Rnw'

###################################################
### code chunk number 1: qualint.Rnw:25-26
###################################################
library("QualInt")


###################################################
### code chunk number 2: qualint.Rnw:33-38
###################################################
ynorm <- rnorm(300)
trtment <- sample(c(0, 1), 300, prob = c(0.4, 0.6), 
                  replace = TRUE)
subgrp <- sample(c(0, 1, 2), 300, prob = c(1/3, 1/3, 1/3), 
                 replace = TRUE)


###################################################
### code chunk number 3: qualint.Rnw:43-44
###################################################
out1 <- qualint(ynorm, trtment, subgrp)


###################################################
### code chunk number 4: qualint.Rnw:49-50
###################################################
out2 <- qualint(ynorm, trtment, subgrp, type = "continuous", plotout = TRUE)


###################################################
### code chunk number 5: qualint.Rnw:55-56
###################################################
out3 <- qualint(ynorm, trtment, subgrp, test = "LRT", type = "continuous")


###################################################
### code chunk number 6: qualint.Rnw:61-62
###################################################
print(out1)


###################################################
### code chunk number 7: qualint.Rnw:67-68
###################################################
print(out3)


###################################################
### code chunk number 8: qualint.Rnw:73-74
###################################################
summary(out1)


###################################################
### code chunk number 9: qualint.Rnw:79-80
###################################################
out1$pvalue


###################################################
### code chunk number 10: qualint.Rnw:91-92
###################################################
coef(out1)


###################################################
### code chunk number 11: qualint.Rnw:98-99
###################################################
plot(out1)


###################################################
### code chunk number 12: qualint.Rnw:103-104
###################################################
ibga(out1)


###################################################
### code chunk number 13: qualint.Rnw:111-114
###################################################
ybin <- sample(c(0, 1), 300, prob = c(0.3, 0.7), replace = TRUE)
out4 <- qualint(ybin, trtment, subgrp, type = "binary", scale = "RD")
print(out4)


###################################################
### code chunk number 14: qualint.Rnw:119-122
###################################################
out5 <- qualint(ybin, trtment, subgrp, type = "binary", 
                scale = "RR", test = "LRT")
print(out5)


###################################################
### code chunk number 15: qualint.Rnw:131-135
###################################################
time <- rpois(300, 200)
censor <- sample(c(0, 1), 300, prob = c(0.7, 0.3), replace = TRUE)
out6 <- qualint(Surv(time, censor), trtment, subgrp, scale = "HR", type = "survival")
print(out6)


