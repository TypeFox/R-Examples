### R code from vignette source 'multinomial-addiction1.Rnw'

###################################################
### code chunk number 1: multinomial-addiction1.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=85)


###################################################
### code chunk number 2: multinomial-addiction1.Rnw:19-22 (eval = FALSE)
###################################################
## library(catdata)
## data(addiction)
## attach(addiction)


###################################################
### code chunk number 3: multinomial-addiction1.Rnw:27-28 (eval = FALSE)
###################################################
## library(nnet)


###################################################
### code chunk number 4: multinomial-addiction1.Rnw:33-35 (eval = FALSE)
###################################################
## ill <- as.factor(ill)
## addiction$ill<-as.factor(addiction$ill)


###################################################
### code chunk number 5: multinomial-addiction1.Rnw:40-42 (eval = FALSE)
###################################################
## multinom0 <- multinom(ill ~ gender + age + university, data=addiction)
## summary(multinom0)


###################################################
### code chunk number 6: multinomial-addiction1.Rnw:46-50 (eval = FALSE)
###################################################
## library(VGAM)
## multivgam0<-vglm(ill ~ gender + age + university, multinomial(refLevel=1), 
##                  data=addiction)
## summary(multivgam0)


###################################################
### code chunk number 7: multinomial-addiction1.Rnw:56-63 (eval = FALSE)
###################################################
## addiction$age2 <- addiction$age^2
## multinom1 <- update(multinom0, . ~ . + age2)
## summary(multinom1)
## 
## multivgam1<-vglm(ill ~ gender + age + university + age2, multinomial(refLevel=1), 
##                  data=addiction)
## summary(multivgam1)


###################################################
### code chunk number 8: multinomial-addiction1.Rnw:69-71 (eval = FALSE)
###################################################
## anova(multinom0,multinom1)
## multinom1$dev - multinom0$dev


###################################################
### code chunk number 9: multinomial-addiction1.Rnw:76-81 (eval = FALSE)
###################################################
## minage <- min(na.omit(age))
## maxage <- max(na.omit(age))
## 
## ageindex <- seq(minage, maxage, 0.1)
## n <- length(ageindex)


###################################################
### code chunk number 10: multinomial-addiction1.Rnw:87-97 (eval = FALSE)
###################################################
## ageindex2 <- ageindex^2
## 
## gender1 <- rep(1, n)
## gender0 <- rep(0, n)
## university1 <- rep(1, n)
## 
## datamale <- as.data.frame(cbind(gender=gender0,age=ageindex,university=
##   university1,age2=ageindex2)) 
## datafemale <- as.data.frame(cbind(gender=gender1,age=ageindex,university=
##   university1,age2=ageindex2))


###################################################
### code chunk number 11: multinomial-addiction1.Rnw:102-104 (eval = FALSE)
###################################################
## probsmale <- predict(multinom1, datamale, type="probs")
## probsfemale <- predict(multinom1, datafemale, type="probs")


###################################################
### code chunk number 12: multinomial-addiction1.Rnw:109-117 (eval = FALSE)
###################################################
## par(cex=1.4, lwd=2)
## 
## plot(ageindex, probsmale[,1], type="l", lty=1, ylim=c(0,1), main=
## "men with university degree", ylab="probabilities")
## lines(ageindex, probsmale[,2], lty="dotted") 
## lines(ageindex, probsmale[,3], lty="dashed") 
## legend("topright", legend=c("Weak-willed", "diseased", "both"), lty=c("solid", 
## "dotted", "dashed"))


###################################################
### code chunk number 13: multinomial-addiction1.Rnw:120-128 (eval = FALSE)
###################################################
## par(cex=1.4, lwd=2)
## 
## plot(ageindex, probsfemale[,1], type="l", lty=1, ylim=c(0,1), main=
##   "women with university degree", ylab="probabilities")
## lines(ageindex, probsfemale[,2], lty="dotted")                                                  
## lines(ageindex, probsfemale[,3], lty="dashed") 
## legend("topright", legend=c("Weak-willed", "diseased", "both"), 
##        lty=c("solid", "dotted", "dashed"))


