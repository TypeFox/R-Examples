
setMethod("plot", signature(x="ur.ers", y="missing"), function(x){
  if(is.null(x@testreg)){
    stop("No plot method for P-test available")}
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,1))
  layout(matrix(c(1, 2, 3, 1, 2, 4), 3 , 2))
  suppressWarnings(plot.ts(diff(x@yd)[-c(1:x@lag)], main="Diagram of fit for test regression", sub=paste("detrending ", x@model, " and ", x@lag, " lagged differences used in test regression",  sep=""), ylab="Actual and fitted values", xlab=""))
  lines(diff(x@yd)[-c(1:x@lag)] - resid(x@testreg), col="seagreen")
  plot.ts(resid(x@testreg), main="Residuals", ylab="", xlab="")
  abline(h=0, col="red")
  acf(resid(x@testreg), main="Autocorrelations of Residuals")
  pacf(resid(x@testreg), main="Partial Autocorrelations of Residuals")
})

setMethod("plot", signature(x="ca.jo", y="missing"), function(x){
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(2,1))
  if(x@P==nrow(x@V)){
    ci <- x@x%*%x@V
  }else{
    ci <- x@x%*%x@V[-(x@P+1),]
  }
  for( i in 1:x@P){
    plot.ts(x@x[,i], main=paste("Time series plot of y", i, sep=""), ylab="")
    plot.ts(ci[,i], main=paste("Cointegration relation of ", i, ". variable", sep=""), ylab="")
    if(interactive()){
      cat("\nType <Return> to continue: ")
      readline()
    }
  }
})

setMethod("plot", signature(x="ur.kpss", y="missing"), function(x){
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,1))
  layout(matrix(c(1, 2, 1, 3), 2 , 2))
  plot.ts(x@res, main=paste("Residuals from test regression of type:", x@type, " with", x@lag, "lags", sep=" "), ylab="residuals", xlab="")
  abline(h=0, col="red")
  acf(x@res, main="Autocorrelations of Residuals")
  pacf(x@res, main="Partial Autocorrelations of Residuals")
})

setMethod("plot", signature(x="ca.po", y="missing"), function(x){
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,1))
  if(x@type=="Pu"){
    layout(matrix(c(1, 2, 1, 3), 2 , 2))
    suppressWarnings(plot.ts(x@res[,1], main="Residuals of CI-regression for y1", sub=paste("detrending:", x@model, sep=" "), ylab="", xlab=""))
    abline(h=0, col="red")
    acf(x@res[,1], main="Autocorrelations of Residuals")
    pacf(x@res[,1], main="Partial Autocorrelations of Residuals")
  }else if(x@type=="Pz"){
    m <- ncol(x@z)
    for( i in 1:m){
      layout(matrix(c(1, 2, 1, 3), 2 , 2))
      suppressWarnings(plot.ts(x@res[,i], main=paste("Residuals of CI-regression with y", i, " as lhs", sep=""), sub=paste("detrending:", x@model, sep=" "), ylab="", xlab=""))
      abline(h=0, col="red")
      acf(x@res[,i], main="Autocorrelations of Residuals")
      pacf(x@res[,i], main="Partial Autocorrelations of Residuals")
      if(interactive()){
        cat("\nType <Return> to continue: ")
        readline()
      }
    }
  }     
})

setMethod("plot", signature(x="ur.pp", y="missing"), function(x){
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,1))
  layout(matrix(c(1, 2, 3, 1, 2, 4), 3 , 2))
  plot.ts(x@y[-1], main=paste("Diagram of fit for model", x@model, sep=" "), ylab="Actual and fitted values", xlab="")
  lines(x@y - x@res, col="seagreen")
  plot.ts(x@res, main="Residuals", ylab="", xlab="")
  abline(h=0, col="red")
  acf(x@res, main="Autocorrelations of Residuals")
  pacf(x@res, main="Partial Autocorrelations of Residuals")
})

setMethod("plot", signature(x="ur.df", y="missing"), function(x){
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,1))
  layout(matrix(c(1, 2, 1, 3), 2 , 2))
  plot.ts(x@res, main="Residuals", ylab="", xlab="")
  abline(h=0, col="red")
  acf(x@res, main="Autocorrelations of Residuals")
  pacf(x@res, main="Partial Autocorrelations of Residuals")
})

setMethod("plot", signature(x="ur.sp", y="missing"), function(x){
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,1))
  layout(matrix(c(1, 2, 3, 1, 2, 4), 3 , 2))
  plot.ts(x@y[-1], main=paste("Diagram of fit for model with polynomial degree of ", x@polynomial, sep="") , ylab="Actual and fitted values", xlab="")
  lines(x@y[-1] - x@res, col="seagreen")
  plot.ts(x@res, main="Residuals", ylab="", xlab="")
  abline(h=0, col="red")
  acf(x@res, main="Autocorrelations of Residuals")
  pacf(x@res, main="Partial Autocorrelations of Residuals")
})

setMethod("plot", signature(x="ur.za", y="missing"), function(x){
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,1))
  yvals <- sort(c(x@cval, x@tstats))
  n <- length(x@y)
  xvals <- pretty(1:n)
  plot.ts(x@tstats, main="Zivot and Andrews Unit Root Test", ylab="t-statistics for lagged endogenous variable", ylim=c(min(yvals), max(yvals)))
  abline(h=x@cval, col=c("red", "blue", "seagreen"))
  if(x@teststat < x@cval[3]){
    abline(v=x@bpoint, col="red", lty=2)}
  mtext(paste("Model type:", x@model, sep=" "), side=1, line=4)
  legend(x=n, y=max(yvals), c("1% c.v.", "2.5% c.v.", "5% c.v."), col=c("red", "blue", "seagreen"), xjust=1, yjust=1, lty=1, horiz=TRUE, cex=0.66, bty="n")
})
  
