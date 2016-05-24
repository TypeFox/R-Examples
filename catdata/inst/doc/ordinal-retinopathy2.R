### R code from vignette source 'ordinal-retinopathy2.Rnw'

###################################################
### code chunk number 1: ordinal-retinopathy2.Rnw:13-15 (eval = FALSE)
###################################################
## rm(list=ls())
## options(width=80)


###################################################
### code chunk number 2: ordinal-retinopathy2.Rnw:18-21 (eval = FALSE)
###################################################
## library(catdata)
## data(retinopathy)
## attach(retinopathy)


###################################################
### code chunk number 3: ordinal-retinopathy2.Rnw:25-26 (eval = FALSE)
###################################################
## library(VGAM)


###################################################
### code chunk number 4: ordinal-retinopathy2.Rnw:32-35 (eval = FALSE)
###################################################
## seqm1 <- vglm(RET ~ SM + DIAB + GH + BP, family = sratio (link="logit", 
## parallel=FALSE))
## deviance(seqm1)


###################################################
### code chunk number 5: ordinal-retinopathy2.Rnw:40-43 (eval = FALSE)
###################################################
## seqm2 <- vglm(RET ~ SM + DIAB + GH + BP, family = sratio (link="logit", 
## parallel=FALSE ~ SM + GH + BP))
## deviance(seqm2)


###################################################
### code chunk number 6: ordinal-retinopathy2.Rnw:47-48 (eval = FALSE)
###################################################
## 1-pchisq(deviance(seqm2)-deviance(seqm1), df=1)


###################################################
### code chunk number 7: ordinal-retinopathy2.Rnw:53-56 (eval = FALSE)
###################################################
## seqm3 <- vglm(RET ~ SM + DIAB + GH + BP, family = sratio (link="logit", 
## parallel=FALSE ~ SM + BP))
## deviance(seqm3)


###################################################
### code chunk number 8: ordinal-retinopathy2.Rnw:60-61 (eval = FALSE)
###################################################
## 1-pchisq(deviance(seqm3)-deviance(seqm2), df=1)


###################################################
### code chunk number 9: ordinal-retinopathy2.Rnw:66-69 (eval = FALSE)
###################################################
## seqm4 <- vglm(RET ~ SM + DIAB + GH + BP, family = sratio (link="logit", 
## parallel=FALSE ~ SM))
## deviance(seqm4)


###################################################
### code chunk number 10: ordinal-retinopathy2.Rnw:73-74 (eval = FALSE)
###################################################
## 1-pchisq(deviance(seqm4)-deviance(seqm3), df=1)


###################################################
### code chunk number 11: ordinal-retinopathy2.Rnw:79-82 (eval = FALSE)
###################################################
## seqm5 <- vglm(RET ~ SM + DIAB + GH + BP, family = sratio (link="logit", 
## parallel=TRUE))
## deviance(seqm5)


###################################################
### code chunk number 12: ordinal-retinopathy2.Rnw:86-87 (eval = FALSE)
###################################################
## 1-pchisq(deviance(seqm5)-deviance(seqm4), df=1)


###################################################
### code chunk number 13: ordinal-retinopathy2.Rnw:92-93 (eval = FALSE)
###################################################
## summary(seqm4)


###################################################
### code chunk number 14: ordinal-retinopathy2.Rnw:103-104 (eval = FALSE)
###################################################
## 1 - pchisq(9.5223^2, df=1)


###################################################
### code chunk number 15: ordinal-retinopathy2.Rnw:108-109 (eval = FALSE)
###################################################
## 1 - pchisq(8.9957^2, df=1)


###################################################
### code chunk number 16: ordinal-retinopathy2.Rnw:113-114 (eval = FALSE)
###################################################
## 1 - pchisq((-1.8646)^2, df=1)


###################################################
### code chunk number 17: ordinal-retinopathy2.Rnw:118-119 (eval = FALSE)
###################################################
## 1 - pchisq(1.5687^2, df=1)


###################################################
### code chunk number 18: ordinal-retinopathy2.Rnw:123-124 (eval = FALSE)
###################################################
## 1 - pchisq((-10.4303)^2, df=1)


###################################################
### code chunk number 19: ordinal-retinopathy2.Rnw:128-129 (eval = FALSE)
###################################################
## 1 - pchisq((-6.3116)^2, df=1)


###################################################
### code chunk number 20: ordinal-retinopathy2.Rnw:133-134 (eval = FALSE)
###################################################
## 1 - pchisq((-5.1037)^2, df=1)


###################################################
### code chunk number 21: ordinal-retinopathy2.Rnw:139-140 (eval = FALSE)
###################################################
## exp(coefficients(seqm4)[3:7])


###################################################
### code chunk number 22: ordinal-retinopathy2.Rnw:143-144 (eval = FALSE)
###################################################
## detach(retinopathy)


###################################################
### code chunk number 23: ordinal-retinopathy2.Rnw:147-148 (eval = FALSE)
###################################################
## detach(package:VGAM)


