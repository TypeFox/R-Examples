### R code from vignette source 'ExactBinarySequentialDesigns.Rnw'

###################################################
### code chunk number 1: ExactBinarySequentialDesigns.Rnw:23-24
###################################################
#library(clinfun)


###################################################
### code chunk number 2: ExactBinarySequentialDesigns.Rnw:43-46
###################################################
library(binseqtest)
b0<-designOBFpower(100, theta0=.4,theta1=.8,power = 0.8, tsalpha = c(0.05, 0.025))
stopTable(b0)


###################################################
### code chunk number 3: ExactBinarySequentialDesigns.Rnw:52-53
###################################################
plot(b0)


###################################################
### code chunk number 4: ExactBinarySequentialDesigns.Rnw:75-76
###################################################
B<-designAb(Nk=c(13,43),a=c(3),theta0=.2,conf.level=.95,alternative="greater")


###################################################
### code chunk number 5: ExactBinarySequentialDesigns.Rnw:83-84
###################################################
plot(B)


###################################################
### code chunk number 6: ExactBinarySequentialDesigns.Rnw:90-91
###################################################
powerBsb(B,theta=.4)


###################################################
### code chunk number 7: ExactBinarySequentialDesigns.Rnw:95-96
###################################################
stopTable(B,S=20,N=43)


###################################################
### code chunk number 8: ExactBinarySequentialDesigns.Rnw:100-102
###################################################
B2<-modify(B,conf.level=.95,alternative="two.sided")
stopTable(B2,S=20,N=43)


###################################################
### code chunk number 9: ExactBinarySequentialDesigns.Rnw:106-108
###################################################
B3<-modify(B,missN=13:18)
stopTable(B3,S=20,N=43)


###################################################
### code chunk number 10: ExactBinarySequentialDesigns.Rnw:113-114
###################################################
plot(B3)


###################################################
### code chunk number 11: ExactBinarySequentialDesigns.Rnw:125-126
###################################################
b<-designOBF(50,k=Inf,theta0=.3)


###################################################
### code chunk number 12: ExactBinarySequentialDesigns.Rnw:133-144
###################################################
b<-designOBF(50,k=Inf,theta0=.3)
par(mfrow=c(2,2))
plot(b,bplottype="NS")
title("bplottype=NS")
plot(b,bplottype="FS")
title("bplottype=FS")
plot(b,bplottype="NE")
title("bplottype=NE")
plot(b,bplottype="NB")
title("bplottype=NB")
par(mfrow=c(1,1))


