plot.dma <-
function(x, ...){
  par(mfrow=c(1,2))
  matplot(x$pmp,main="Posterior model probabilities",col=1:ncol(x$pmp),
          ylab="Posterior model probability",xlab="Index",type='l')
  ltmp<-"Model 1"
  if(ncol(x$pmp)>1){
  for(i in 2:ncol(x$pmp)){ltmp<-c(ltmp,paste("Model",i))}}
  legend(x=round(nrow(x$pmp)/2,0),y=round((max(x$pmp/2)),1),col=c(1:3),legend=ltmp,
         lty=c(1:3))
  plot(x$yhat.ma,type='l',main="Model averaged fitted value", ylab="Fitted value")
}

