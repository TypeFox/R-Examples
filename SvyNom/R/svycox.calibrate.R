svycox.calibrate <-
function(.nom,.timept=.nom$pred.at,.ngroup=5){
.loc=max(which(.nom$preds[[1]]$time<=.timept))
pred.timept=rep(NA,.nom$svy.cox$n)
for(i in 1:length(pred.timept)) pred.timept[i]=.nom$preds[[i]]$surv[.loc]
pred.timept.grp=cut(pred.timept,quantile(pred.timept,seq(0,1,1/.ngroup)),labels=1:.ngroup)
.predicted=tapply(pred.timept,pred.timept.grp,median)
.observed=matrix(NA,nrow=.ngroup,ncol=3)
colnames(.observed)=c("Observed","Lower 95%","Upper 95%")
for(i in 1:.ngroup){
  .km1=svykm(as.formula(paste(names(.nom$svy.cox$model)[1],"~","1")), 
			design=subset(.nom$design,pred.timept.grp==i), se=TRUE)
  .km1.timept=.km1[[2]][which(.km1[[1]]>.timept)[1]-1]
  .varlog1.timept=.km1[[3]][which(.km1[[1]]>.timept)[1]-1]
  .ll1.timept=exp(log(.km1.timept)-1.96*sqrt(.varlog1.timept))
  .ul1.timept=exp(log(.km1.timept)+1.96*sqrt(.varlog1.timept))
  .observed[i,]=c(.km1.timept,.ll1.timept,.ul1.timept)
  }
plot(.predicted,.observed[,1],xlim=0:1,ylim=0:1,xlab="Predicted",ylab="Observed",pch=16,type="b")
arrows(x0=.predicted,y0=.observed[,2],y1=.observed[,3],angle=90,code=3,length=0.05,lwd=1)
abline(0,1,lty=2)
return(cbind(Predicted=.predicted,.observed))
}

