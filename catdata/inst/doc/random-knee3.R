### R code from vignette source 'random-knee3.Rnw'

###################################################
### code chunk number 1: random-knee3.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: random-knee3.Rnw:19-21 (eval = FALSE)
###################################################
## library(catdata)
## data(knee)


###################################################
### code chunk number 3: random-knee3.Rnw:28-36 (eval = FALSE)
###################################################
## knee <- reshape(knee, direction="long", varying=list(5:8), v.names="R",
##                 timevar="Time")
## 
## knee$RD <- rep(0, length(knee$R))
## knee$RD[knee$R>2] <- 1
## 
## knee$Age <- knee$Age - 30
## knee$Age2<-knee$Age^2


###################################################
### code chunk number 4: random-knee3.Rnw:41-43 (eval = FALSE)
###################################################
## knee$Th <- as.factor(knee$Th)
## knee$Sex <- as.factor(knee$Sex)


###################################################
### code chunk number 5: random-knee3.Rnw:48-49 (eval = FALSE)
###################################################
## library(flexmix)


###################################################
### code chunk number 6: random-knee3.Rnw:55-63 (eval = FALSE)
###################################################
## kneeflex2 <-stepFlexmix(cbind(RD,1-RD) ~ 1 | N,	model = FLXMRglmfix(family = 
##   "binomial", fixed= ~ Th + Sex + Age + Age2), k = 2, nrep = 5, data = knee)
## 
## kneeflex3 <-stepFlexmix(cbind(RD,1-RD) ~ 1 | N,	model = FLXMRglmfix(family = 
##   "binomial", fixed= ~ Th + Sex + Age + Age2), k = 3, nrep = 5, data = knee)
## 
## kneeflex4 <-stepFlexmix(cbind(RD,1-RD) ~ 1 | N,	model = FLXMRglmfix(family = 
##   "binomial", fixed= ~ Th + Sex + Age + Age2), k = 4, nrep = 5, data = knee)


###################################################
### code chunk number 7: random-knee3.Rnw:69-72 (eval = FALSE)
###################################################
## summary(kneeflex2)@BIC
## summary(kneeflex3)@BIC
## summary(kneeflex4)@BIC


###################################################
### code chunk number 8: random-knee3.Rnw:77-78 (eval = FALSE)
###################################################
## summary(refit(kneeflex2))


###################################################
### code chunk number 9: random-knee3.Rnw:81-82 (eval = FALSE)
###################################################
## detach(package:flexmix)


