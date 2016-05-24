### R code from vignette source 'dlnmOverview.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(continue="  ")


###################################################
### code chunk number 2: load
###################################################
library(dlnm)


###################################################
### code chunk number 3: changelog (eval = FALSE)
###################################################
## file.show(system.file("ChangeLog", package="dlnm"))


###################################################
### code chunk number 4: poly
###################################################
dlnm:::poly(1:5,degree=3)


###################################################
### code chunk number 5: strata
###################################################
dlnm:::strata(1:5,breaks=c(2,4))[,]


###################################################
### code chunk number 6: thr
###################################################
dlnm:::thr(1:5,thr.value=3,side="d")[,]


###################################################
### code chunk number 7: onebasis1
###################################################
onebasis(1:5,fun="poly",degree=3)


###################################################
### code chunk number 8: onebasis2
###################################################
mylog <- function(x) log(x)
onebasis(1:5,"mylog")


###################################################
### code chunk number 9: histories
###################################################
hist <- matrix(sample(1:12),3,dimnames=list(paste("sub",1:3,sep=""),
  paste("lag",2:5,sep="")))
hist


###################################################
### code chunk number 10: crossbasis1
###################################################
crossbasis(hist,lag=c(2,5),argvar=list(fun="poly",degree=2),
  arglag=list(fun="strata",breaks=4))[,]


###################################################
### code chunk number 11: crossbasis
###################################################
cb <- crossbasis(chicagoNMMAPS$temp,lag=30,argvar=list("thr",thr.value=c(10,20)),
  arglag=list(knots=c(1,4,12)))
summary(cb)


###################################################
### code chunk number 12: regression
###################################################
model <- lm(cvd~cb,chicagoNMMAPS)
pred <- crosspred(cb,model,at=-20:30)


###################################################
### code chunk number 13: prediction1
###################################################
c(pred$matfit["-10","lag5"],pred$matlow["-10","lag5"],pred$mathigh["-10","lag5"])
pred$allfit["25"]


###################################################
### code chunk number 14: prediction2
###################################################
histpred <- t(c(rep(30,10),rep(22,21)))
crosspred(cb,model,at=histpred)$allfit


###################################################
### code chunk number 15: plot3dcontournoeval (eval = FALSE)
###################################################
## plot(pred,ptype="3d",main="3D plot",xlab="Temperature",zlab="CVD excess count",
##   theta=200, ltheta=180)
## plot(pred,ptype="contour",key.title=title("CVD"),
##   plot.title=title("Contour plot",xlab="Temperature",ylab="Lag"))


###################################################
### code chunk number 16: plot3d
###################################################
plot(pred,ptype="3d",main="3D plot",xlab="Temperature",zlab="CVD excess count",
  theta=200, ltheta=180)


###################################################
### code chunk number 17: plotcontour
###################################################
plot(pred,ptype="contour",key.title=title("CVD"),
  plot.title=title("Contour plot",xlab="Temperature",ylab="Lag"))


###################################################
### code chunk number 18: plotoverall
###################################################
plot(pred,"overall",xlab="Temperature",ylab="CVD counts",ylim=c(-5,25))


###################################################
### code chunk number 19: plotslicenoeval (eval = FALSE)
###################################################
## plot(pred,"slices",lag=5,xlab="Temperature",ylab="CVD count",
##   col="blue",ci="lines",ci.arg=list(lty=5))
## plot(pred,"slices",var=25,xlab="Lag",ylab="CVD count",type="p",pch=19,ci="bars")


###################################################
### code chunk number 20: plotslicelag
###################################################
plot(pred,"slices",lag=5,xlab="Temperature",ylab="CVD count",
  col="blue",ci="lines",ci.arg=list(lty=5))


###################################################
### code chunk number 21: plotslicevar
###################################################
plot(pred,"slices",var=25,xlab="Lag",ylab="CVD count",type="p",pch=19,ci="bars")


###################################################
### code chunk number 22: internal1 (eval = FALSE)
###################################################
## help(getcoef)


###################################################
### code chunk number 23: internal1 (eval = FALSE)
###################################################
## dlnm:::fci
## getAnywhere(fci)


###################################################
### code chunk number 24: changelog (eval = FALSE)
###################################################
## file.show(system.file("ChangeLog", package="dlnm"))
## file.show(system.file("Changesince151", package="dlnm"))
## file.show(system.file("Changesince200", package="dlnm"))
## file.show(system.file("Changesince220", package="dlnm"))


