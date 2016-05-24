### R code from vignette source 'irtguide.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: irtguide.Rnw:68-69
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: irtguide.Rnw:200-232
###################################################
library(kequate)
set.seed(7)
akX <- runif(15, 0.5, 2)
bkX <- rnorm(15)
ckX <- runif(15, 0.1, 0.2)
akY <- runif(15, 0.5, 2)
bkY <- rnorm(15)
ckY <- runif(15, 0.1, 0.2)
akA <- runif(15, 0.5, 2)
bkA <- rnorm(15)
ckA <- runif(15, 0.1, 0.2)

dataP <- matrix(0, nrow=1000, ncol=30)
dataQ <- matrix(0, nrow=1000, ncol=30)
data3plP <- matrix(0, nrow=1000, ncol=30)
data3plQ <- matrix(0, nrow=1000, ncol=30)

for(i in 1:1000){
ability <- rnorm(1)
dataP[i,1:15] <- (1/(1+exp(-akX*(ability-bkX)))) > runif(15)
dataP[i,16:30] <- (1/(1+exp(-akA*(ability-bkA)))) > runif(15)
data3plP[i,1:15] <- (ckX+(1-ckX)/(1+exp(-akX*(ability-bkX)))) > runif(15)
data3plP[i,16:30] <- (ckA+(1-ckA)/(1+exp(-akA*(ability-bkA)))) > runif(15)
}

for(i in 1:1000){
ability <- rnorm(1, mean=0.5)
dataQ[i,1:15] <- (1/(1+exp(-akY*(ability -bkY)))) > runif(15)
dataQ[i,16:30] <- (1/(1+exp(-akA*(ability -bkA)))) > runif(15)
data3plQ[i,1:15] <- (ckY+(1-ckY)/(1+exp(-akY*(ability-bkY)))) > runif(15)
data3plQ[i,16:30] <- (ckA+(1-ckA)/(1+exp(-akA*(ability-bkA)))) > runif(15)
}


###################################################
### code chunk number 3: irtguide.Rnw:236-237
###################################################
eq2pl <- irtose("CE", dataP, dataQ, 0:15, 0:15, 0:15)


###################################################
### code chunk number 4: irtguide.Rnw:240-241
###################################################
summary(eq2pl)


###################################################
### code chunk number 5: irtguide.Rnw:246-247
###################################################
irtobjects <- eq2pl@irt


###################################################
### code chunk number 6: irtguide.Rnw:250-252
###################################################
sim2plP <- irtobjects$ltmP
sim2plQ <- irtobjects$ltmQ


###################################################
### code chunk number 7: irtguide.Rnw:256-257
###################################################
load("irtguide.RData")


###################################################
### code chunk number 8: irtguide.Rnw:260-262
###################################################
eq3pl <- irtose("CE", sim3plP, sim3plQ, 0:15, 0:15, 0:15, model="3pl")
summary(eq3pl)


###################################################
### code chunk number 9: eq3plplot
###################################################
plot(eq3pl)


###################################################
### code chunk number 10: eq3plplot1
###################################################
plot(eq3pl)


###################################################
### code chunk number 11: irtguide.Rnw:285-288
###################################################
eq2plLOW <- irtose("CE", sim2plP, sim2plQ, 0:15, 0:15, 0:15, qpoints=-1)
eq2plAVG <- irtose("CE", sim2plP, sim2plQ, 0:15, 0:15, 0:15, qpoints=0)
eq2plHIGH <- irtose("CE", sim2plP, sim2plQ, 0:15, 0:15, 0:15, qpoints=1)


###################################################
### code chunk number 12: leplot
###################################################
plot(0:15, getEq(eq2plLOW), ylim=c(-1, 15), pch=1, ylab="", xlab="")
par(new=TRUE)
plot(0:15, getEq(eq2plAVG), ylim=c(-1, 15), pch=2, ylab="", xlab="")
par(new=TRUE)
plot(0:15, getEq(eq2plHIGH), ylim=c(-1, 15), pch=3, ylab="Equated value", xlab="Score value")
legend("topleft", inset=.1, title="Ability level:", c("-1", "0", "1"), pch=c(1, 2, 3))


###################################################
### code chunk number 13: leplot1
###################################################
plot(0:15, getEq(eq2plLOW), ylim=c(-1, 15), pch=1, ylab="", xlab="")
par(new=TRUE)
plot(0:15, getEq(eq2plAVG), ylim=c(-1, 15), pch=2, ylab="", xlab="")
par(new=TRUE)
plot(0:15, getEq(eq2plHIGH), ylim=c(-1, 15), pch=3, ylab="Equated value", xlab="Score value")
legend("topleft", inset=.1, title="Ability level:", c("-1", "0", "1"), pch=c(1, 2, 3))


