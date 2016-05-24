
## ----mandn,warning=FALSE,message=FALSE,fig.width=15,results='hide'-------
library(knitr)
library(hysteresis)
par(mfrow=c(3,3),mai=c(0,0.2,0.2,0),ann=FALSE,xaxt="n",yaxt="n",oma=c(0,0,3,0))

for (i in c(1,3,15)){
  for (j in c(1,3,15)){
    obj<-mloop(m=i,n=j,n.points=100,period=99)
    plot(floop(obj$x,obj$y,m=i,n=j,period=99),xlim=c(-0.8,0.8),ylim=c(-0.8,0.8))
    if (i==1) title(paste("n=",j,sep=""))
    if (j==1) title(ylab=paste("m=",i,sep=""),line=0,cex.sub=2)
  }
}
title("Hysteresis Loops for Odd Values of m and n",outer=TRUE)


## ----moremandn,warning=FALSE, fig.width=15-------------------------------

par(mfrow=c(3,3),mai=c(0,0.2,0.2,0),ann=FALSE,xaxt="n",yaxt="n",oma=c(0,0,3,0))

for (i in c(1,3,15)){
  for (j in c(2,4,16)){
    obj<-mloop(m=i,n=j,n.points=100,period=99)
    plot(floop(obj$x,obj$y,m=i,n=j,period=99),xlim=c(-0.8,0.8),ylim=c(-0.8,0.8))
    if (i==1) title(paste("n=",j,sep=""))
    if (j==2) title(ylab=paste("m=",i,sep=""),line=0,cex.sub=2)
  }
}
title("Hysteresis Loops for Odd Values of m and Even Values of n",outer=TRUE)


## ------------------------------------------------------------------------
obj<-mloop(cx=0,cy=0,n.points=100,period=99)
obj2<-mloop(cx=1.5,cy=0,n.points=100,period=99)
obj3<-mloop(cx=0,cy=1.5,n.points=100,period=99)
plot(obj$x,obj$y,type="l",xlim=c(-2,3),ylim=c(-2,3),xlab="x",ylab="y",col="#6600CC",main="Centroid Given by cx and cy")
points(0,0,pch=19,col="#6600CC")
text(x=0,y=0.15,"(cx=0,cy=0)",col="#6600CC")
lines(obj2$x,obj2$y,col="#00FF66")
points(1.5,0,pch=19,col="#00FF66")
text(x=1.5,y=0.15,"(cx=1.5,cy=0)",col="#00FF66")
lines(obj3$x,obj3$y,col="#FF6600")
points(0,1.5,pch=19,col="#FF6600")
text(x=0,y=1.65,"(cx=0,cy=1.5)",col="#FF6600")


## ----bx,fig.show='asis',fig.width=5--------------------------------------
for (i in c(1,2,4)){
  obj<-mloop(b.x=i,n.points=100,period=99)
  plot(obj$x,obj$y,xlim=c(-5,10),ylim=c(-1.4,1.4),type="l",main=paste("b.x=",i,sep=""),xlab="x",ylab="y")
  points(i,0.8,pch=19)
  legend(i,1,legend=c("Saturation Point","x=cx+b.x","y=cy+b.y"),bty="n")
}


## ----by,fig.show='asis',fig.width=5--------------------------------------
for (i in c(0.8,2,4)){
  obj<-mloop(b.y=i,n.points=100,period=99)
  plot(obj$x,obj$y,xlim=c(-1,2),ylim=c(-5,5),type="l",main=paste("b.y=",i,sep=""),xlab="x",ylab="y")
  points(0.6,i,pch=19)
  legend(0.6,i,legend=c("Saturation Point","x=cx+b.x","y=cy+b.y"),bty="n")
}


## ----retention,fig.show='asis',fig.width=5-------------------------------
for (i in c(1,2,4)){
  obj<-mloop(retention=i,n.points=100,period=99)
  plot(obj$x,obj$y,xlim=c(-1,1),ylim=c(-5,5),type="l",main=paste("retention=",i,sep=""),xlab="x",ylab="y")
  segments(0,0,0,i)
  text(0.3,0.5,"Retention")
}


## ------------------------------------------------------------------------
obj<-mloop(retention=0.5,n.points=100,period=99)
  plot(obj$x,obj$y,type="l",xlab="x",ylab="y",main="Starting Points for Different Values of phase.angle",xlim=c(-0.6,0.8))
for (i in c(0,90,180,260)){
  obj2<-mloop(phase.angle=i,retention=0.5,n.points=1,period=99)
  points(obj2$x,obj2$y,pch=19,col="gold",cex=2)
  points(obj2$x,obj2$y,col="gold",cex=4)
  text(obj2$x+.08,obj2$y,round(i,2))
}


## ------------------------------------------------------------------------
set.seed(24)
ellipse1 <- mel(method=2,retention=0.4,b.x=0.6,b.y=0.8,cx=0,cy=0,sd.x=0.1,sd.y=0.1,phase.angle=0,period=24,n.points=24)
#The function **mel** can be used as an alternative to **mloop** for simulating ellipses, and it is useful because it offers four different ellipse parameterizations.
model <- fel(ellipse1$x,ellipse1$y,method="harmonic2",period=24,times="equal")
#period=24 and times="equal" are used to say that 24 equally spaced points make up an ellipse.
model


## ------------------------------------------------------------------------
model$Estimates


## ------------------------------------------------------------------------
summary(model,N=10000,studentize=TRUE)


## ------------------------------------------------------------------------

plot(model,main="2-step Simple Harmonic Regression Ellipse Example")


## ------------------------------------------------------------------------
modeldirect <- fel(ellipse1$x,ellipse1$y,method="direct",period=24,times="equal")
summodel<-summary(modeldirect,N=10000,studentize=TRUE)
summodel
plot(modeldirect,main="Direct Ellipse Example")


## ------------------------------------------------------------------------
summodel$values


## ----tests,fig.width=10,fig.height=10,warning=FALSE,message=FALSE--------
set.seed(13)
par(mfrow=c(2,2))
halfellipse <- mel(method=2,cx=20,cy=25,retention=1.2,b.x=14,b.y=0.8,sd.x=1,sd.y=0.2,period=24,n.points=16,phase.angle=pi/2)
halftrueellipse <- mel(method=2,cx=20,cy=25,retention=1.2,b.x=14,b.y=0.8,phase.angle=0,period=99,n.points=100)
harmodel<-fel(halfellipse$x,halfellipse$y,method="harmonic2",period=24,times="equal")
dirmodel<-fel(halfellipse$x,halfellipse$y,method="direct",period=24,times="equal")
lmmodel<-fel(halfellipse$x,halfellipse$y,method="lm",period=24,times="equal")
nlsmodel<-fel(halfellipse$x,halfellipse$y,method="nls",period=24,times="equal",control=c(n.iter=500))
plot(harmodel,main="Harmonic2 Model",xlim=c(5,34),ylim=c(23.4,26.9))
lines(halftrueellipse$x,halftrueellipse$y,col="red")
plot(dirmodel,main="Direct Model",xlim=c(5,34),ylim=c(23.4,26.9))
lines(halftrueellipse$x,halftrueellipse$y,col="red")
plot(lmmodel,main="Linear Model",xlim=c(5,34),ylim=c(23.4,26.9))
lines(halftrueellipse$x,halftrueellipse$y,col="red")
plot(nlsmodel,main="Non-Linear Model",xlim=c(5,34),ylim=c(23.4,26.9))
lines(halftrueellipse$x,halftrueellipse$y,col="red")


## ----geometric,warning=FALSE---------------------------------------------
set.seed(101)
ellip <- mel(rote.deg=45,semi.major=5,semi.minor=3,n.points=13,sd.x=0.4,sd.y=0.4)
true.ellip <- mel(rote.deg=45,semi.major=5,semi.minor=3,n.points=100,period=100)
ellip.geometric <- fel(ellip$x,ellip$y,method="geometric")
ellip.geometric$values
plot(ellip.geometric,main="Geometric Model")
lines(true.ellip$x,true.ellip$y,col="red")


## ----boottest,fig.width=10,fig.height=10---------------------------------
par(mfrow=c(2,2))
harsummodel<-summary(harmodel,N=1000,studentize=TRUE)
dirsummodel<-summary(dirmodel,N=1000,studentize=TRUE)
lmsummodel<-summary(lmmodel,N=1000,studentize=TRUE)
nlssummodel<-summary(nlsmodel,N=1000,studentize=TRUE)
plot(harsummodel,main="Bootstrapped Harmonic2 Model",xlim=c(5,34),ylim=c(23.4,26.9))
lines(halftrueellipse$x,halftrueellipse$y,col="red")
plot(dirsummodel,main="Bootstrapped Direct Model",xlim=c(5,34),ylim=c(23.4,26.9))
lines(halftrueellipse$x,halftrueellipse$y,col="red")
plot(lmsummodel,main="Bootstrapped Lm Model",xlim=c(5,34),ylim=c(23.4,26.9))
lines(halftrueellipse$x,halftrueellipse$y,col="red")
plot(nlssummodel,main="Bootstrapped Nls Model",xlim=c(5,34),ylim=c(23.4,26.9))
lines(halftrueellipse$x,halftrueellipse$y,col="red")


## ----multiple,warning=FALSE----------------------------------------------
data(EllipseData)
models <- fel(EllipseData$X,EllipseData$Y,method="harmonic2",subjects=EllipseData$subjects,subset=EllipseData$repeated==1)
models
summodels<-summary(models)
summodels
plot(summodels,main="Fitting Multiple Ellipses Simultaneously")


## ----mid,eval=FALSE------------------------------------------------------
## ## write.table(models$Estimates,"file_name.txt") and
## ## write.table(summodels$Boot.Estimates,"file_name.txt")


## ------------------------------------------------------------------------
loop <- mloop(n=5, m=3,sd.x=0.02,sd.y=0.02)
fitloop <- floop(loop$x,loop$y,n=5, m=3,period=24,times="equal")
fitloop$Estimates
plot(fitloop,main="Fitted Hysteresis Loop")
summary(fitloop)


