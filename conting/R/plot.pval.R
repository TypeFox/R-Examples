plot.pval <-
function(x,...){

plot(x$Tobs,x$Tpred,xlab="T_obs",ylab="T_pred",type="n",...)
points(x$Tobs[x$Tpred>x$Tobs],x$Tpred[x$Tpred>x$Tobs],col=8,pch=16)
points(x$Tobs[x$Tpred<x$Tobs],x$Tpred[x$Tpred<x$Tobs],pch=16)
abline(a=0,b=1,lty=2,lwd=2)
legend(x=min(x$Tobs),y=max(x$Tpred),legend=c("T_pred>T_obs","T_pred<T_obs"),col=c(8,1),pch=c(16,16))

}
