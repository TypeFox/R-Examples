library(conicfit)

x <- matrix(c(1, 2, 5, 7, 9, 6, 3, 8, 7, 6, 8, 7, 5, 7, 2, 4),8,2,byrow=FALSE)

xl<-c(-2,10)
yl<-c(-2,12)
plot(x[,1],x[,2],col='red',xlim=xl,ylim=yl,xlab='',ylab='');par(new=TRUE)

centroid<-apply(x,2,mean)
x<-x-matrix(rep(centroid,dim(x)[1]),dim(x)[1],2,byrow=T)

b<-fitbookstein(x)
xyDirect<-calculateEllipse(b[["z"]][1]+centroid[1], b[["z"]][2]+centroid[2],b[["a"]], b[["b"]],  180/pi*b[["alpha"]])
plot(xyDirect[,1],xyDirect[,2],xlim=xl,ylim=yl,type='l',col='blue',xlab='',ylab='');par(new=TRUE)

b<-fitggk(x)
xyDirect<-calculateEllipse(b[["z"]][1]+centroid[1], b[["z"]][2]+centroid[2],b[["a"]], b[["b"]],  180/pi*b[["alpha"]])
plot(xyDirect[,1],xyDirect[,2],xlim=xl,ylim=yl,type='l',col='green',xlab='',ylab='');par(new=TRUE)

ellipDirect <- EllipseDirectFit(x)

estCircle <- estimateInitialGuessCircle(x)
elliLMG <- fit.ellipseLMG(x,cbind(c(estCircle,estCircle[3],0)),0.1)$ParG
xyLMG <-calculateEllipse(elliLMG[1]+centroid[1], elliLMG[2]+centroid[2], elliLMG[3], elliLMG[4], 180/pi*elliLMG[5])
plot(xyLMG[,1],xyLMG[,2],xlim=xl,ylim=yl,type='l',col='orange',xlab='',ylab='');par(new=TRUE)

elliLMA <- fit.conicLMA(x,ellipDirect,0.1)$ParA
elliLMAG <- AtoG(elliLMA)$ParG
xyLMA<-calculateEllipse(elliLMAG[1]+centroid[1], elliLMAG[2]+centroid[2], elliLMAG[3], elliLMAG[4], 180/pi*elliLMAG[5])
plot(xyLMA[,1],xyLMA[,2],xlim=xl,ylim=yl,type='l',col='purple',xlab='',ylab='',lty=2);par(new=TRUE)

legend('topleft', c('Data points','Algebraic linear - Bookstein','Algebraic linear - Trace','NLS LMA','NLS LMG'), text.col=c('red','blue','green','orange','purple'))