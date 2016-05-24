plotSEPmvr <-
function(mvrdcvobj,optcomp,y,X,method="simpls",complete=TRUE, ...)
{
# Generate plot showing SEP values for
# Repeated Double Cross Validation 
#

dat=list(y=y,X=as.matrix(X))
SEP.dcv=apply(drop(mvrdcvobj$pred-drop(y)),c(3,2),sd) # sd not from matrix!
matplot(t(SEP.dcv),type="l",lty=1,col=gray(0.6),
  xlab="Number of components",ylab="SEP",cex.lab=1.2,...)

#require(pls)
if (complete){
  ncomp <- ncol(SEP.dcv)
}
else {
  ncomp <- optcomp
}
mvr.cv=mvr(y~X,ncomp=ncomp,data=dat,method=method,validation="CV")
mvr.res=(drop(mvr.cv$vali$pred)-drop(y))
mvr.SEP=apply(mvr.res,2,sd)
lines(mvr.SEP,col=1,lwd=2)

abline(v=optcomp,lty=2)
abline(h=mvrdcvobj$SEPopt,lty=2)

list(SEPdcv=SEP.dcv,SEPcv=mvr.SEP)
}

