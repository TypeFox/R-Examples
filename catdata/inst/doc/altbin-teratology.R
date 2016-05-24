### R code from vignette source 'altbin-teratology.Rnw'

###################################################
### code chunk number 1: altbin-teratology.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=60)


###################################################
### code chunk number 2: altbin-teratology.Rnw:20-23 (eval = FALSE)
###################################################
## library(catdata)
## data(teratology)
## data(teratology2)


###################################################
### code chunk number 3: altbin-teratology.Rnw:28-29 (eval = FALSE)
###################################################
## attach(teratology)


###################################################
### code chunk number 4: altbin-teratology.Rnw:33-35 (eval = FALSE)
###################################################
## mLogit <- glm(cbind(D,L) ~ as.factor(Grp), family=binomial())
## summary(mLogit)


###################################################
### code chunk number 5: altbin-teratology.Rnw:41-43 (eval = FALSE)
###################################################
## mQuasi <- glm(cbind(D,L) ~ as.factor(Grp), family=quasibinomial(link="logit"))
## summary(mQuasi)


###################################################
### code chunk number 6: altbin-teratology.Rnw:49-50 (eval = FALSE)
###################################################
## library(gee)


###################################################
### code chunk number 7: altbin-teratology.Rnw:55-57 (eval = FALSE)
###################################################
## detach(teratology)
## attach(teratology2)


###################################################
### code chunk number 8: altbin-teratology.Rnw:64-66 (eval = FALSE)
###################################################
## mGee <- gee(y ~ as.factor(Grp), id=Rat, family=binomial)
## summary(mGee) 


###################################################
### code chunk number 9: altbin-teratology.Rnw:71-72 (eval = FALSE)
###################################################
## library(VGAM)


###################################################
### code chunk number 10: altbin-teratology.Rnw:75-77 (eval = FALSE)
###################################################
## detach(teratology2)
## attach(teratology)


###################################################
### code chunk number 11: altbin-teratology.Rnw:83-84 (eval = FALSE)
###################################################
## N <- D + L


###################################################
### code chunk number 12: altbin-teratology.Rnw:89-91 (eval = FALSE)
###################################################
## mBetaBin <- vglm(cbind(D,L) ~ as.factor(Grp), family=betabinomial, subset=N>1)
## summary(mBetaBin)


###################################################
### code chunk number 13: altbin-teratology.Rnw:96-98 (eval = FALSE)
###################################################
## detach(teratology)
## attach(teratology2)


###################################################
### code chunk number 14: altbin-teratology.Rnw:104-106 (eval = FALSE)
###################################################
## mMixPql<- glmmPQL(y ~ as.factor(Grp), random=~1 | Rat, family=binomial)
## summary(mMixPql)


###################################################
### code chunk number 15: altbin-teratology.Rnw:111-112 (eval = FALSE)
###################################################
## library(glmmML)


###################################################
### code chunk number 16: altbin-teratology.Rnw:118-121 (eval = FALSE)
###################################################
## mGaussH <- glmmML(y ~ as.factor(Grp), cluster=Rat, method = "ghq", n.points = 14,
##                   boot = 0) 
## summary(mGaussH)


###################################################
### code chunk number 17: altbin-teratology.Rnw:126-128 (eval = FALSE)
###################################################
## detach(teratology2)
## attach(teratology)


###################################################
### code chunk number 18: altbin-teratology.Rnw:133-134 (eval = FALSE)
###################################################
## library(flexmix)


###################################################
### code chunk number 19: altbin-teratology.Rnw:142-145 (eval = FALSE)
###################################################
## detach(package:nlme)
## detach(package:VGAM)
## library(stats4)


###################################################
### code chunk number 20: altbin-teratology.Rnw:148-152 (eval = FALSE)
###################################################
## mDiscmix <-stepFlexmix(cbind(D,L) ~ 1, k = 2, nrep=5,
##                  model = FLXMRglmfix(family = "binomial",fixed =~as.factor(Grp)))
## summary(mDiscmix)
## parameters(mDiscmix)


