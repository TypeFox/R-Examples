plotarimapred <-
function(ts.cont, fit.arima, xlim, range.percent=0.2, xreg=NULL, ylab=NULL, xlab=NULL, main=NULL){
  if(!is.null(xreg)) ts.pred <- forecast(fit.arima,length(ts.cont),xreg=xreg)
  else ts.pred <- forecast(fit.arima,length(ts.cont))
  
  mn <- min(ts.pred$lower,ts.cont)
  mx <- max(ts.pred$upper,ts.cont)
  
  yrange.min <- mn-range.percent*abs(mn)
  yrange.max <- mx+range.percent*abs(mx)
  yrange <- range(yrange.min,yrange.max)
  
  plot(ts.pred,xlim=xlim,ylim=yrange, ylab=ylab, xlab=xlab, main=main)
  lines(ts.cont,xlim=xlim, lty=5, lwd=2, col="black")
}
