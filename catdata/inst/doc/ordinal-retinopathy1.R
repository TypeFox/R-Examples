### R code from vignette source 'ordinal-retinopathy1.Rnw'

###################################################
### code chunk number 1: ordinal-retinopathy1.Rnw:12-14 (eval = FALSE)
###################################################
## options(width=80)
## rm(list=ls(all=TRUE))


###################################################
### code chunk number 2: ordinal-retinopathy1.Rnw:17-20 (eval = FALSE)
###################################################
## library(catdata)
## data(retinopathy)
## attach(retinopathy)


###################################################
### code chunk number 3: ordinal-retinopathy1.Rnw:28-31 (eval = FALSE)
###################################################
## library(VGAM)
## RET <- as.ordered(RET)
## SM <- as.factor(SM)


###################################################
### code chunk number 4: ordinal-retinopathy1.Rnw:37-38 (eval = FALSE)
###################################################
## pom <- vglm(RET ~ SM + DIAB + GH + BP, family = cumulative (parallel=TRUE))


###################################################
### code chunk number 5: ordinal-retinopathy1.Rnw:41-42 (eval = FALSE)
###################################################
## ppom <- vglm(RET ~ SM + DIAB + GH + BP, family = cumulative (parallel=FALSE))


###################################################
### code chunk number 6: ordinal-retinopathy1.Rnw:47-48 (eval = FALSE)
###################################################
## deviance(pom)


###################################################
### code chunk number 7: ordinal-retinopathy1.Rnw:51-52 (eval = FALSE)
###################################################
## deviance(ppom)


###################################################
### code chunk number 8: ordinal-retinopathy1.Rnw:57-58 (eval = FALSE)
###################################################
## 1 - pchisq(deviance(pom) - deviance(ppom), df=4)


###################################################
### code chunk number 9: ordinal-retinopathy1.Rnw:65-66 (eval = FALSE)
###################################################
## summary(pom)


###################################################
### code chunk number 10: ordinal-retinopathy1.Rnw:71-72 (eval = FALSE)
###################################################
## summary(ppom)


###################################################
### code chunk number 11: ordinal-retinopathy1.Rnw:83-86 (eval = FALSE)
###################################################
## ppom2 <- vglm (RET ~ SM + DIAB + GH + BP,
## family = cumulative (parallel = FALSE ~ SM + DIAB + GH))
## deviance(ppom2)


###################################################
### code chunk number 12: ordinal-retinopathy1.Rnw:89-90 (eval = FALSE)
###################################################
## 1-pchisq(deviance(ppom2)-deviance(ppom), df=1)


###################################################
### code chunk number 13: ordinal-retinopathy1.Rnw:95-98 (eval = FALSE)
###################################################
## ppom3 <- vglm (RET ~ SM + DIAB + GH + BP, 
## family = cumulative (parallel = FALSE ~ SM + DIAB))
## deviance(ppom3)


###################################################
### code chunk number 14: ordinal-retinopathy1.Rnw:101-102 (eval = FALSE)
###################################################
## 1-pchisq(deviance(ppom3)-deviance(ppom2), df=1)


###################################################
### code chunk number 15: ordinal-retinopathy1.Rnw:106-109 (eval = FALSE)
###################################################
## ppom4 <- vglm (RET ~ SM + DIAB + GH + BP, 
## family = cumulative (parallel = FALSE ~ SM))
## deviance(ppom4)


###################################################
### code chunk number 16: ordinal-retinopathy1.Rnw:112-113 (eval = FALSE)
###################################################
## 1-pchisq(deviance(ppom4)-deviance(ppom3), df=1)


###################################################
### code chunk number 17: ordinal-retinopathy1.Rnw:117-118 (eval = FALSE)
###################################################
## 1-pchisq(deviance(pom)-deviance(ppom4), df=1)


###################################################
### code chunk number 18: ordinal-retinopathy1.Rnw:121-122 (eval = FALSE)
###################################################
## detach(retinopathy)


###################################################
### code chunk number 19: ordinal-retinopathy1.Rnw:124-125 (eval = FALSE)
###################################################
## detach(package:VGAM)


