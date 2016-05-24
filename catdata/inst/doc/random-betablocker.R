### R code from vignette source 'random-betablocker.Rnw'

###################################################
### code chunk number 1: random-betablocker.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: random-betablocker.Rnw:19-21 (eval = FALSE)
###################################################
## library(flexmix)
## data(betablocker)


###################################################
### code chunk number 3: random-betablocker.Rnw:24-25 (eval = FALSE)
###################################################
## betablocker$Treatment <- as.factor(betablocker$Treatment)


###################################################
### code chunk number 4: random-betablocker.Rnw:30-33 (eval = FALSE)
###################################################
## GlmT <- glm(cbind(Deaths, Total	- Deaths) ~ Treatment, family = "binomial", 
## data = betablocker)
## summary(GlmT)


###################################################
### code chunk number 5: random-betablocker.Rnw:39-42 (eval = FALSE)
###################################################
## GlmTC <- glm(cbind(Deaths, Total - Deaths) ~ Treatment + as.factor(Center), 
##              family 	= "binomial", data = betablocker)
## summary(GlmTC)


###################################################
### code chunk number 6: random-betablocker.Rnw:47-48 (eval = FALSE)
###################################################
## library(glmmML)


###################################################
### code chunk number 7: random-betablocker.Rnw:53-57 (eval = FALSE)
###################################################
## MixedGH4 <- glmmML(cbind(Deaths, Total - Deaths) ~ Treatment, cluster=Center, 
##                    method = c("ghq"), n.points = 4, boot = 0, data=betablocker)
## 
## summary(MixedGH4)


###################################################
### code chunk number 8: random-betablocker.Rnw:62-66 (eval = FALSE)
###################################################
## MixedGH20 <- glmmML(cbind(Deaths, Total - Deaths) ~ Treatment, cluster=Center, 
##                     method = c("ghq"), n.points = 20, boot = 0, data=betablocker)
## 
## summary(MixedGH20)


###################################################
### code chunk number 9: random-betablocker.Rnw:69-70 (eval = FALSE)
###################################################
## set.seed(5)


###################################################
### code chunk number 10: random-betablocker.Rnw:76-77 (eval = FALSE)
###################################################
## detach(package:glmmML)


###################################################
### code chunk number 11: random-betablocker.Rnw:80-83 (eval = FALSE)
###################################################
## MixFix3 <-stepFlexmix(cbind(Deaths, Total - Deaths) ~ 1 | Center,	model = 
##   FLXMRglmfix(family = "binomial", fixed = ~ Treatment), k = 3, nrep = 5, 
##                       data = betablocker)


###################################################
### code chunk number 12: random-betablocker.Rnw:88-89 (eval = FALSE)
###################################################
## MixFix3


###################################################
### code chunk number 13: random-betablocker.Rnw:94-95 (eval = FALSE)
###################################################
## parameters(MixFix3)


###################################################
### code chunk number 14: random-betablocker.Rnw:100-101 (eval = FALSE)
###################################################
## library(flexmix)


###################################################
### code chunk number 15: random-betablocker.Rnw:104-106 (eval = FALSE)
###################################################
## summary(MixFix3)
## summary(refit(MixFix3))


###################################################
### code chunk number 16: random-betablocker.Rnw:109-110 (eval = FALSE)
###################################################
## set.seed(5)


###################################################
### code chunk number 17: random-betablocker.Rnw:115-118 (eval = FALSE)
###################################################
## MixFix4 <-stepFlexmix(cbind(Deaths, Total - Deaths) ~ 1 | Center, model = 
##   FLXMRglmfix(family = "binomial", fixed = ~ Treatment), k = 4, nrep = 5, 
##                       data = betablocker)


###################################################
### code chunk number 18: random-betablocker.Rnw:121-122 (eval = FALSE)
###################################################
## MixFix4


###################################################
### code chunk number 19: random-betablocker.Rnw:125-126 (eval = FALSE)
###################################################
## parameters(MixFix4)


###################################################
### code chunk number 20: random-betablocker.Rnw:129-131 (eval = FALSE)
###################################################
## summary(MixFix4)
## summary(refit(MixFix4))


###################################################
### code chunk number 21: random-betablocker.Rnw:133-134 (eval = FALSE)
###################################################
## detach(package:flexmix)


