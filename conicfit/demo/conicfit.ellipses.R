# start the RNG with a seed for reproducibility
set.seed(0)
# 50 points from an ellipse at c(0,0) with axis (200, 100), angle 45 degrees
a<-calculateEllipse(0,0,200,100,45,50)
plot(a[,1],a[,2],xlim=c(-250,250),ylim=c(-250,250),type='l');par(new=TRUE)
plot(0,0,xlim=c(-250,250),ylim=c(-250,250),pch=4);par(new=TRUERUE) # plot the center with an X
# 50 points from an ellipse at c(0,0) with axis (200, 100), angle 45 degrees, positioned by uniform random distribution, noise=normal random distribution with sd=50
xy<-calculateEllipse(0,0,200,100,45,50, randomDist=TRUE,noiseFun=function(x) (x+rnorm(1,mean=0,sd=50)))
plot(xy[,1],xy[,2],xlim=c(-250,250),ylim=c(-250,250),col='magenta');par(new=TRUE)

ellipTaubin <- EllipseFitByTaubin(xy)
ellipTaubinG <- AtoG(ellipTaubin)$ParG
xyTaubin<-calculateEllipse(ellipTaubinG[1], ellipTaubinG[2], ellipTaubinG[3], ellipTaubinG[4], 180/pi*ellipTaubinG[5])
plot(xyTaubin[,1],xyTaubin[,2],xlim=c(-250,250),ylim=c(-250,250),type='l',col='red');par(new=TRUE)

ellipDirect <- EllipseDirectFit(xy)
ellipDirectG <- AtoG(ellipDirect)$ParG
xyDirect<-calculateEllipse(ellipDirectG[1], ellipDirectG[2], ellipDirectG[3], ellipDirectG[4], 180/pi*ellipDirectG[5])
plot(xyDirect[,1],xyDirect[,2],xlim=c(-250,250),ylim=c(-250,250),type='l',col='cyan');par(new=TRUE)

estCircle <- estimateInitialGuessCircle(xy)

elliLMG <- fit.ellipseLMG(xy,cbind(c(estCircle,estCircle[3],0)),0.1)$ParG
xyLMG <-calculateEllipse(elliLMG[1], elliLMG[2], elliLMG[3], elliLMG[4], 180/pi*elliLMG[5])
plot(xyLMG[,1],xyLMG[,2],xlim=c(-250,250),ylim=c(-250,250),type='l',col='orange');par(new=TRUE)

elliLMA <- fit.conicLMA(xy,ellipDirect,0.1)$ParA
elliLMAG <- AtoG(elliLMA)$ParG
xyLMA<-calculateEllipse(elliLMAG[1], elliLMAG[2], elliLMAG[3], elliLMAG[4], 180/pi*elliLMAG[5])
plot(xyLMA[,1],xyLMA[,2],xlim=c(-250,250),ylim=c(-250,250),type='l',col='purple');par(new=TRUE)

