plot.BchronRSLRun <-
function(x, xlab='Age (cal BP)',ylab='Depth (m)', ...) {
  
  age.low = apply(x$BchronologyRun$thetaPredict,2,'quantile',probs=0.025)
  age.med = apply(x$BchronologyRun$thetaPredict,2,'quantile',probs=0.5)
  age.high = apply(x$BchronologyRun$thetaPredict,2,'quantile',probs=0.975)
  
  RSL.low = x$RSLmean-2*x$RSLsd
  RSL.high = x$RSLmean+2*x$RSLsd
  
  graphics::plot(age.med,x$RSLmean,type='n',xlim=rev(range(c(age.low,age.high))),ylim=range(c(RSL.low,RSL.high)),xlab=xlab,ylab=ylab,...)
  
  for(i in 1:length(x$RSLmean)) {
    agescale = (age.high[i]-age.low[i])/4
    rslscale = x$RSLsd[i]
    graphics::lines(ellipse::ellipse(0,scale=c(agescale,rslscale),centre=c(age.med[i],x$RSLmean[i])),col=grDevices::rgb(0,0,1,0.4))
  }
  
  xgrid = seq(max(age.low),min(age.high),length=100)/1000
  
  pred.lines = matrix(NA,ncol=length(xgrid),nrow=nrow(x$samples))
  degmat = matrix(rep(0:(x$degree),length(xgrid)*(x$degree+1)),nrow=length(xgrid),ncol=x$degree+1,byrow=TRUE)
  X.pred = matrix(rep(xgrid-x$const,x$degree+1),ncol=x$degree+1)
  X.pred = (X.pred^degmat)
  
  for(i in 1:nrow(pred.lines)) {
    pred.lines[i,] = X.pred%*%matrix(x$samples[i,],ncol=1,nrow=x$degree+1)
  }
  
  pred.med = apply(pred.lines,2,'quantile',probs=0.5)
  pred.low = apply(pred.lines,2,'quantile',probs=0.025)
  pred.high = apply(pred.lines,2,'quantile',probs=0.975)
  graphics::lines(xgrid*1000,pred.med,lwd=2)
  graphics::lines(xgrid*1000,pred.low,lwd=2,lty=2)
  graphics::lines(xgrid*1000,pred.high,lwd=2,lty=2)
  
}
