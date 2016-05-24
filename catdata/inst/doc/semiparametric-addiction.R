### R code from vignette source 'semiparametric-addiction.Rnw'

###################################################
### code chunk number 1: semiparametric-addiction.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: semiparametric-addiction.Rnw:19-22 (eval = FALSE)
###################################################
## library(catdata)
## data(addiction)
## attach(addiction)


###################################################
### code chunk number 3: semiparametric-addiction.Rnw:27-28 (eval = FALSE)
###################################################
## library(mgcv)


###################################################
### code chunk number 4: semiparametric-addiction.Rnw:34-38 (eval = FALSE)
###################################################
## minage <- min(na.omit(age))
## maxage <- max(na.omit(age))
## ageindex <- seq(minage, maxage, 0.1)
## n <- length(ageindex)


###################################################
### code chunk number 5: semiparametric-addiction.Rnw:41-49 (eval = FALSE)
###################################################
## gender1 <- rep(1, n)
## gender0 <- rep(0, n)
## university1 <- rep(1, n)
## university0 <- rep(0, n)
## 
## datamale <- as.data.frame(cbind(gender=gender0,age=ageindex,university=university1))
## datafemale <- as.data.frame(cbind(gender=gender1,age=ageindex,
##                                   university=university1))


###################################################
### code chunk number 6: semiparametric-addiction.Rnw:54-57 (eval = FALSE)
###################################################
## ill01 <-ill
## ill01[ill==1] <- 0
## ill01[ill==2] <- 1


###################################################
### code chunk number 7: semiparametric-addiction.Rnw:63-66 (eval = FALSE)
###################################################
## gam1 <- gam(as.factor(ill01) ~ s(age) + gender + university, family=binomial())
## gam2 <- gam(as.factor(ill) ~ s(age) + gender + university, family=binomial(), 
##             data=addiction[addiction$ill!=2,])


###################################################
### code chunk number 8: semiparametric-addiction.Rnw:71-73 (eval = FALSE)
###################################################
## summary(gam1)
## summary(gam2)


###################################################
### code chunk number 9: semiparametric-addiction.Rnw:78-79 (eval = FALSE)
###################################################
## plot(gam1)


###################################################
### code chunk number 10: semiparametric-addiction.Rnw:82-83 (eval = FALSE)
###################################################
## plot(gam2)


###################################################
### code chunk number 11: semiparametric-addiction.Rnw:89-91 (eval = FALSE)
###################################################
## predmale1 <- predict(gam1, newdata=datamale, type="response")
## predmale2 <- predict(gam2, newdata=datamale, type="response")


###################################################
### code chunk number 12: semiparametric-addiction.Rnw:94-96 (eval = FALSE)
###################################################
## predfemale1 <- predict(gam1, newdata=datafemale, type="response")
## predfemale2 <- predict(gam2, newdata=datafemale, type="response")


###################################################
### code chunk number 13: semiparametric-addiction.Rnw:101-104 (eval = FALSE)
###################################################
## p2 <- predmale1
## p1 <- predmale2 * (1-predmale1)
## p0 <- (1-predmale2) * (1-predmale1)


###################################################
### code chunk number 14: semiparametric-addiction.Rnw:107-110 (eval = FALSE)
###################################################
## pf2 <- predfemale1
## pf1 <- predfemale2 * (1-predfemale1)
## pf0 <- (1-predfemale2) * (1-predfemale1)


###################################################
### code chunk number 15: semiparametric-addiction.Rnw:115-130 (eval = FALSE)
###################################################
## par(mfrow=c(1,2), cex=1.8)
## plot(ageindex, p0, type="l", lty=1, ylim=c(0,1), main="men with university degree",
##      ylab="probabilities")
## lines(ageindex, p1, lty="dotted")
## lines(ageindex, p2, lty="dashed")
## legend("topright", legend=c("Weak-willed", "diseased", "both"), 
##        lty=c("solid", "dotted", "dashed"))
## 
## 
## plot(ageindex, pf0, type="l", lty=1, ylim=c(0,1), 
##      main="women with university degree", ylab="probabilities")
## lines(ageindex, pf1, lty="dotted")
## lines(ageindex, pf2, lty="dashed")
## legend("topright", legend=c("Weak-willed", "diseased", "both"),
##        lty=c("solid", "dotted", "dashed"))


###################################################
### code chunk number 16: semiparametric-addiction.Rnw:135-139 (eval = FALSE)
###################################################
## gam3 <- gam(as.factor(ill)~ s(age) + gender + university, 
##             data=addiction[addiction$ill!=2,], family=binomial())
## gam4 <- gam(as.factor(ill)~ s(age) + gender + university, 
##             data=addiction[addiction$ill!=1,], family=binomial())


###################################################
### code chunk number 17: semiparametric-addiction.Rnw:143-145 (eval = FALSE)
###################################################
## summary(gam3)
## summary(gam4)


