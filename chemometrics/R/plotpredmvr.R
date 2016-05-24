plotpredmvr <-
function(mvrdcvobj,optcomp,y,X,method="simpls",...)
{
# Generate plot showing predicted values for
# Repeated Double Cross Validation 
#

dat=list(y=y,X=as.matrix(X))
#require(pls)
mvr.cv=mvr(y~X,ncomp=optcomp,data=dat,method=method,validation="CV")

par(mfrow=c(1,2))
ylimits=range(mvrdcvobj$pred[,,optcomp,])
plot(y,mvr.cv$vali$pred[,,optcomp],xlab="Measured y",ylab="Predicted y",
        cex.lab=1.2,cex=0.7,pch=3,col=1,ylim=ylimits,...)
title("Prediction from CV")
abline(c(0,1),lty=1)
plot(y,apply(mvrdcvobj$pred[,,optcomp,],1,mean),xlab="Measured y",
        ylab="Predicted y",cex.lab=1.2,type="n",ylim=ylimits,...)
for (i in 1:ncol(mvrdcvobj$pred[,,optcomp,])){
  points(y,mvrdcvobj$pred[,,optcomp,i],pch=3,cex=0.7,col=gray(0.6))
}
points(y,apply(mvrdcvobj$pred[,,optcomp,],1,mean),pch=3,cex=0.7,col=1)
title("Prediction from Repeated Double-CV")
abline(c(0,1),lty=1)

invisible()
}

