plotresmvr <-
function(mvrdcvobj,optcomp,y,X,method="simpls",...)
{
# Generate plot showing residuals for
# Repeated Double Cross Validation
#

dat=list(y=y,X=as.matrix(X))
#require(pls)
mvr.cv=mvr(y~X,ncomp=optcomp,data=dat,method=method,validation="CV")

par(mfrow=c(1,2))
predcv=mvr.cv$vali$pred[,,optcomp]
preddcvall=mvrdcvobj$pred[,,optcomp,]
preddcv=apply(preddcvall,1,mean)
ylimits=max(abs(preddcvall-drop(y)))
ylimits=sort(c(-ylimits,ylimits))
plot(predcv,predcv-y,xlab="Predicted y",ylab="Residuals",cex.lab=1.2,
        cex=0.7,pch=3,col=1,ylim=ylimits,...)
title("Results from CV")
abline(h=0,lty=1)
plot(preddcv,preddcv-y,xlab="Predicted y",ylab="Residuals",cex.lab=1.2,
        cex=0.7,pch=3,col=gray(0.6),type="n",ylim=ylimits,...)
for (i in 1:ncol(preddcvall)){
  points(preddcv,preddcvall[,i]-y,cex=0.7,pch=3,col=gray(0.6))
}
points(preddcv,preddcv-y,cex=0.7,pch=3,col=1)
title("Results from Repeated Double-CV")
abline(h=0,lty=1)

invisible()
}

