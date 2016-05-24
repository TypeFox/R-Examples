### R code from vignette source 'dbmss.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Declarations
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("dbmss")


###################################################
### code chunk number 2: SimpleWmppp
###################################################
Pattern <- wmppp(data.frame(X=runif(100), Y=runif(100)))
summary(Pattern)


###################################################
### code chunk number 3: CSBIGS
###################################################
load("CSBIGS.Rdata")
Category <- cut(Emergencies$M, quantile(Emergencies$M, c(0, 0.9, 1)),
   labels = c("Other", "Biggest"), include.lowest = TRUE)
X <- wmppp(data.frame(X=Emergencies$X, Y=Emergencies$Y, PointType=Category),
   win=Region)
X$window$units <- c("meter","meters")


###################################################
### code chunk number 4: Toulouse
###################################################
load("CSBIGS.Rdata")
Category <- cut(Emergencies$M, quantile(Emergencies$M, c(0, 0.9, 1)),
   labels = c("Other", "Biggest"), include.lowest = TRUE)
X <- wmppp(data.frame(X=Emergencies$X, Y=Emergencies$Y, PointType=Category),
   win=Region)
X$window$units <- c("meter","meters")
X2 <- split(X)
marks(X2$Other) <- rep(1, X2$Other$n)
marks(X2$Biggest) <- rep(1, X2$Biggest$n)
par(mfrow=c(1,2), mar=c(0,0,0,0))  
plot(X2$Other, main="", maxsize=1, legend=FALSE)
text(514300,  1826800, "a")
plot(X2$Biggest, main="", maxsize=1, legend=FALSE)
text(514300,  1826800, "b")
par(mfrow=c(1,1))


###################################################
### code chunk number 5: KdCode (eval = FALSE)
###################################################
## load("CSBIGS.Rdata")
## Category <- cut(Emergencies$M, quantile(Emergencies$M, c(0, 0.9, 1)),
##    labels = c("Other", "Biggest"), include.lowest = TRUE)
## X <- wmppp(data.frame(X=Emergencies$X, Y=Emergencies$Y, PointType=Category),
##    win=Region)
## X$window$units <- c("meter","meters")
## KdE <- KdEnvelope(X, r=seq(0, 10000, 100), NumberOfSimulations=1000,
##    ReferenceType="Biggest", Global=TRUE)
## plot(KdE, main="")


###################################################
### code chunk number 6: KdFig (eval = FALSE)
###################################################
## load("CSBIGS.Rdata")
## Category <- cut(Emergencies$M, quantile(Emergencies$M, c(0, 0.9, 1)),
##    labels = c("Other", "Biggest"), include.lowest = TRUE)
## X <- wmppp(data.frame(X=Emergencies$X, Y=Emergencies$Y, PointType=Category),
##    win=Region)
## X$window$units <- c("meter","meters")
## KdE <- KdEnvelope(X, r=seq(0, 10000, 100), NumberOfSimulations=1000,
##    ReferenceType="Biggest", Global=TRUE)
## plot(KdE, main="")


###################################################
### code chunk number 7: P16Code (eval = FALSE)
###################################################
## data("paracou16")
## plot(paracou16, which.marks="PointWeight", main="", legend=FALSE)


###################################################
### code chunk number 8: P16Fig
###################################################
par(mar=c(0,1,0,4))
data("paracou16")
plot(paracou16, which.marks="PointWeight", main="", legend=FALSE)


###################################################
### code chunk number 9: MCode (eval = FALSE)
###################################################
## Envelope <- MEnvelope(paracou16, r = seq(0, 30, 2), NumberOfSimulations 
##    = 1000, Alpha = 0.05, ReferenceType = "V. Americana", NeighborType 
##    = "Q. Rosea", SimulationType = "RandomLabeling", Global = TRUE)
## plot(Envelope, main="", ylim=c(0, 20))


###################################################
### code chunk number 10: MFig
###################################################
Envelope <- MEnvelope(paracou16, r = seq(0, 30, 2), NumberOfSimulations 
   = 1000, Alpha = 0.05, ReferenceType = "V. Americana", NeighborType 
   = "Q. Rosea", SimulationType = "RandomLabeling", Global = TRUE)
plot(Envelope, main="", ylim=c(0, 20))


###################################################
### code chunk number 11: MEnvelope
###################################################
GoFtest(Envelope)


###################################################
### code chunk number 12: Ktest
###################################################
data("paracou16")
RectWindow <- owin(c(300, 400), c(0, 150))
X <- paracou16[RectWindow]
Ktest(X, seq(5, 50, 5))


