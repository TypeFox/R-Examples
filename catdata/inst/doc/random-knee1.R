### R code from vignette source 'random-knee1.Rnw'

###################################################
### code chunk number 1: random-knee1.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: random-knee1.Rnw:20-22 (eval = FALSE)
###################################################
## library(catdata)
## data(knee)


###################################################
### code chunk number 3: random-knee1.Rnw:31-35 (eval = FALSE)
###################################################
## knee <- reshape(knee, direction="long", varying=list(5:8), v.names="R",timevar="Time")
## 
## knee$RD <- rep(0, length(knee$R))
## knee$RD[knee$R>2] <- 1


###################################################
### code chunk number 4: random-knee1.Rnw:40-42 (eval = FALSE)
###################################################
## knee$Age <- knee$Age - 30
## knee$Age2 <- knee$Age^2  


###################################################
### code chunk number 5: random-knee1.Rnw:47-48 (eval = FALSE)
###################################################
## knee <- knee[knee$Time!=1,]


###################################################
### code chunk number 6: random-knee1.Rnw:53-54 (eval = FALSE)
###################################################
## library(glmmML)


###################################################
### code chunk number 7: random-knee1.Rnw:59-62 (eval = FALSE)
###################################################
## kneeGHQ <- glmmML(RD ~ as.factor(Th) + as.factor(Sex) + Age + Age2, data=knee, 
## family=binomial(), method="ghq", n.points=20, cluster=id)
## summary(kneeGHQ)


###################################################
### code chunk number 8: random-knee1.Rnw:67-70 (eval = FALSE)
###################################################
## kneePQL <- glmmPQL(RD ~ as.factor(Th) + as.factor(Sex) + Age + Age2, data=knee, 
## random = ~ 1|id, family=binomial())
## summary(kneePQL)


###################################################
### code chunk number 9: random-knee1.Rnw:75-76 (eval = FALSE)
###################################################
## library(gee)


###################################################
### code chunk number 10: random-knee1.Rnw:82-86 (eval = FALSE)
###################################################
## knee <- knee[order(knee$id),]
## kneeGEE <- gee(RD ~ as.factor(Th) + as.factor(Sex) + Age + Age2, data=knee, 
## family=binomial(), id=id, corstr="exchangeable")
## summary(kneeGEE)


###################################################
### code chunk number 11: random-knee1.Rnw:89-90 (eval = FALSE)
###################################################
## detach(package:gee)


