### R code from vignette source 'ordinal-knee2.Rnw'

###################################################
### code chunk number 1: ordinal-knee2.Rnw:12-14 (eval = FALSE)
###################################################
## options(width=80)
## rm(list=ls(all=TRUE))


###################################################
### code chunk number 2: ordinal-knee2.Rnw:18-21 (eval = FALSE)
###################################################
## library(catdata)
## data(knee)
## attach(knee)


###################################################
### code chunk number 3: ordinal-knee2.Rnw:27-29 (eval = FALSE)
###################################################
## Age <- Age - 30
## Age2 <- Age^2


###################################################
### code chunk number 4: ordinal-knee2.Rnw:35-37 (eval = FALSE)
###################################################
## k <- length(unique(R4))
## R4rev <- k + 1 - R4


###################################################
### code chunk number 5: ordinal-knee2.Rnw:42-46 (eval = FALSE)
###################################################
## R4 <- as.ordered(R4)
## R4rev <- as.ordered(R4rev)
## Th <- as.factor(Th)
## Sex <- as.factor(Sex)


###################################################
### code chunk number 6: ordinal-knee2.Rnw:50-51 (eval = FALSE)
###################################################
## library(VGAM)


###################################################
### code chunk number 7: ordinal-knee2.Rnw:60-62 (eval = FALSE)
###################################################
## clogit <- vglm(R4 ~ Th + Sex + Age +Age2, family = cumulative (link="logit", 
## parallel=TRUE))


###################################################
### code chunk number 8: ordinal-knee2.Rnw:66-68 (eval = FALSE)
###################################################
## cprobit <- vglm(R4 ~ Th + Sex + Age +Age2, family = cumulative (link="probit",
## parallel=TRUE))


###################################################
### code chunk number 9: ordinal-knee2.Rnw:72-74 (eval = FALSE)
###################################################
## cgumbel <- vglm(R4rev ~ Th + Sex + Age +Age2, family = cumulative(link="cloglog", 
## parallel=TRUE))


###################################################
### code chunk number 10: ordinal-knee2.Rnw:78-80 (eval = FALSE)
###################################################
## cgompertz <- vglm(R4 ~ Th + Sex + Age +Age2, family = cumulative(link="cloglog", 
## parallel=TRUE))


###################################################
### code chunk number 11: ordinal-knee2.Rnw:86-90 (eval = FALSE)
###################################################
## deviance(clogit)
## deviance(cprobit)
## deviance(cgumbel)
## deviance(cgompertz)


###################################################
### code chunk number 12: ordinal-knee2.Rnw:98-100 (eval = FALSE)
###################################################
## slogit <- vglm(R4 ~ Th + Sex + Age +Age2, family = sratio (link="logit", 
## parallel=TRUE))


###################################################
### code chunk number 13: ordinal-knee2.Rnw:104-106 (eval = FALSE)
###################################################
## sprobit <- vglm(R4 ~ Th + Sex + Age +Age2, family = sratio (link="probit", 
## parallel=TRUE))


###################################################
### code chunk number 14: ordinal-knee2.Rnw:110-112 (eval = FALSE)
###################################################
## sgumbel <- vglm(R4rev ~ Th + Sex + Age +Age2, family = sratio (link="cloglog", 
## parallel=TRUE))


###################################################
### code chunk number 15: ordinal-knee2.Rnw:116-118 (eval = FALSE)
###################################################
## sgompertz <- vglm(R4 ~ Th + Sex + Age +Age2, family = sratio (link="cloglog",
## parallel=TRUE))


###################################################
### code chunk number 16: ordinal-knee2.Rnw:122-126 (eval = FALSE)
###################################################
## deviance(slogit)
## deviance(sprobit)
## deviance(sgumbel)
## deviance(sgompertz)


###################################################
### code chunk number 17: ordinal-knee2.Rnw:129-130 (eval = FALSE)
###################################################
## detach(package:VGAM)


