### R code from vignette source 'Stationary.rnw'

###################################################
### code chunk number 1: Stationary.rnw:15-19
###################################################
library(TSTutorial)
data(AirBcn)
data(Victimes)
data(Turismes)


###################################################
### code chunk number 2: Stationary.rnw:41-43
###################################################
series=diff(diff(log(AirBcn),12))
ts.plot(series,main="Stationary")


###################################################
### code chunk number 3: Stationary.rnw:54-56
###################################################
series=diff(diff(AirBcn,12))/100
ts.plot(series,main="Nonconstant variance")


###################################################
### code chunk number 4: Stationary.rnw:64-66
###################################################
series=1:300
ts.plot(series,main="Nonconstant mean")


###################################################
### code chunk number 5: Stationary.rnw:74-76
###################################################
series=Turismes
ts.plot(series,main="Seasonal component")


###################################################
### code chunk number 6: Stationary.rnw:84-86
###################################################
series=(1:300)^2
ts.plot(series,main="Nonconstant mean and variance")


###################################################
### code chunk number 7: Stationary.rnw:94-96
###################################################
series=Victimes
ts.plot(series,main="Nonconstant variance and seasonal component")


###################################################
### code chunk number 8: Stationary.rnw:104-106
###################################################
series=log(AirBcn)
ts.plot(series,main="Nonconstant mean and seasonal component")


###################################################
### code chunk number 9: Stationary.rnw:114-116
###################################################
series=AirBcn
ts.plot(series,main="Nonconstant mean and variance, and seasonal comp.")


