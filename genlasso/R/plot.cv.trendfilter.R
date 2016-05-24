plot.cv.trendfilter <- function(x, legendpos="top", xlab, ylab, ...) {
  if (!any(class(x)=="cv.trendfilter")) {
    stop("Passed x must be of class \"cv.trendfilter\".")
  }
  
  if (x$mode=="lambda") {
    xvals = x$lambda
    if (missing(xlab)) xlab = expression(lambda)
    imin = which.min(abs(xvals-x$lambda.min))
    i1se = which.min(abs(xvals-x$lambda.1se))
    log = "x"
  }
  else {
    xvals = x$df
    if (missing(xlab)) xlab = "df"
    imin = which.min(abs(xvals-x$df.min))
    i1se = which.min(abs(xvals-x$df.1se))
    log = ""
  }
  if (missing(ylab)) ylab = "CV error"
  
  plot(xvals,x$err,type="l",xlab=xlab,ylab="CV error",log=log,
       ylim=range(c(x$err+x$se,x$err-x$se)),...)
  lines(xvals,x$err+x$se,lty=3)
  lines(xvals,x$err-x$se,lty=3)
  points(xvals[imin],x$err[imin],pch=19,col="blue")
  points(xvals[i1se],x$err[i1se],pch=19,col="red")
  legend(legendpos, legend=c("Minimal CV error","One standard error rule"),
         col=c("blue","red"), pch=19)
}

