arimapred <-
function(timeseries,timeseries.cont=NULL, n.ahead=NULL, na.action=na.omit, xreg=NULL, newxreg=NULL, se.fit=FALSE, plot=FALSE, range.p=0.2,ylab=NULL,xlab=NULL,main=NULL){
  if(is.null(timeseries))    stop("timeseries is required")
  if(is.null(timeseries.cont) & is.null(n.ahead)) stop("the number of values to be predicted is unknown")
  
  ts <- ts(na.action(timeseries))
  
  N <- NULL
  if(!is.null(timeseries.cont)) N <- length(na.action(timeseries.cont))
  if(!is.null(n.ahead)) N <- n.ahead
  
  nobs <- length(ts)
  i <- nobs + 1
  f <- nobs +  N
  
  reg <- cbind(1:nobs,xreg)
  
  fit <- auto.arima(ts,xreg=ts(reg,start=1))
  
  newreg <- cbind(i:f,newxreg)
  
  pred <- predict(fit,n.ahead=N, newxreg=ts(newreg,start=i),se.fit=se.fit)
  
  if(!is.null(timeseries.cont) & plot)
  {
    ts.cont <- ts(na.action(timeseries.cont),start=i)
    plotarimapred(ts.cont, fit, xlim=c(i,f), range.p, xreg=newreg, ylab=ylab, xlab=xlab, main=main)
  }
  
  return (pred)
}
