### R code from vignette source 'multinomial-party2.Rnw'

###################################################
### code chunk number 1: multinomial-party2.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: multinomial-party2.Rnw:20-32 (eval = FALSE)
###################################################
## partypref <- matrix(data=c(114, 10, 53, 224,134,9,42,226,114,8,23,174,339,30,13,
## 414,42,5,44,161,88,10,60,171,90,8,31,168, 413,23,14,375), nrow=8, byrow=TRUE)
## 
## partydat<-data.frame(
## party=c(rep("CDU",sum(partypref[,1])),rep("SPD",sum(partypref[,4])),
## rep("The Liberals",sum(partypref[,2])),rep("The Greens",sum(partypref[,3]))),
## sex=c(rep(0,sum(partypref[1:4,1])),rep(1,sum(partypref[5:8,1])),
## rep(0,sum(partypref[1:4,4])),rep(1,sum(partypref[5:8,4])),
## rep(0,sum(partypref[1:4,2])),rep(1,sum(partypref[5:8,2])),
## rep(0,sum(partypref[1:4,3])),rep(1,sum(partypref[5:8,3]))),
## age=c(rep(c(1:4,1:4), partypref[,1]),rep(c(1:4,1:4), partypref[,4]),
## rep(c(1:4,1:4), partypref[,2]),rep(c(1:4,1:4), partypref[,3])))


###################################################
### code chunk number 3: multinomial-party2.Rnw:40-47 (eval = FALSE)
###################################################
## x1 <- partypref/rowSums(partypref)
## 
## x1 <- x1[c(1,4,5,8),]
## 
## x1 <- cbind(x1[,2],x1[,3],x1[,1],x1[,4])
## 
## x1 <- x1*6


###################################################
### code chunk number 4: multinomial-party2.Rnw:52-54 (eval = FALSE)
###################################################
## loc1 <- matrix(data=c(0,0,15,0,0,15,15,15),ncol=2,byrow=T)
## library(grDevices)


###################################################
### code chunk number 5: multinomial-party2.Rnw:59-62 (eval = FALSE)
###################################################
## stars(x1, scale=FALSE,key.loc = c(25,6),len=2,cex=1.2,lwd=1.2,
## xlim=c(-10,30),ylim=c(-5,10),key.labels=c("Liberals","Green","CDU","SPD"), 
## location=loc1, labels = c("male young","male old","female young","female old"))


###################################################
### code chunk number 6: multinomial-party2.Rnw:68-75 (eval = FALSE)
###################################################
## partydat$age <- as.factor(partydat$age)
## 
## library(nnet)
## 
## partymult <- multinom(party ~ sex + age , data=partydat)
## 
## summary(partymult)


###################################################
### code chunk number 7: multinomial-party2.Rnw:81-86 (eval = FALSE)
###################################################
## x2 <- matrix(data=c(rep(1,5),as.matrix(t(exp(coefficients(partymult))))),nrow=5)
## 
## x2 <- cbind(x2[,4],x2[,3],x2[,1],x2[,2])
## 
## loc2 <- matrix(data=c(0,9,9,9,0,0,9,0,18,0),ncol=2,byrow=T)


###################################################
### code chunk number 8: multinomial-party2.Rnw:91-94 (eval = FALSE)
###################################################
## stars(x2, scale=FALSE,key.loc = c(18,9),len=2,cex=1.2,lwd=1.2,
## xlim=c(-2,23),ylim=c(-3,10),key.labels=c("Liberals","Green","CDU","SPD"), 
## location=loc2, labels = c("reference","gender","age 2","age 3","age 4"))


