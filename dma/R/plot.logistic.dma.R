plot.logistic.dma <-
function(x, ...){
  par(mfrow=c(1,2))
  matplot(t(x$pmp),main="Posterior model probabilities",col=1:nrow(x$pmp),
          ylab="Posterior model probability",xlab="Index",type='l')
  ltmp<-"Model 1"
  if(ncol(x$pmp)>1){
  for(i in 2:nrow(x$pmp)){ltmp<-c(ltmp,paste("Model",i))}}
  legend(x=0,y=round((max(x$pmp/2)),1),col=c(1:nrow(x$pmp)),legend=ltmp,
         lty=c(1:nrow(x$pmp)))
  plot(x$yhatdma,type='l',main="Model averaged fitted value", ylab="Fitted value")
  #points(lowess(ldma.test$yhatbma~c(1:ncol(x$pmp))),col='red',lwd=2,type='l')
}

