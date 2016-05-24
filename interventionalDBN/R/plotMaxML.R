plotMaxML <-
function(output,xlab="Prior strength",ylab="Marginal likelihood",col.max="red",lty.max=3,lwd.max=1,...) {
  log.ml.prod <- apply(output$marginal.likelihood,2,sum)
  plot(output$inputs$priorStrength,log.ml.prod,t="l",ylim=c(min(log.ml.prod),max(log.ml.prod)),xlab=xlab,ylab=ylab,...)
  segments(output$ebPriorStrength,min(log.ml.prod),output$ebPriorStrength,max(log.ml.prod),col=col.max,lty=lty.max,lwd=lwd.max)  
}
