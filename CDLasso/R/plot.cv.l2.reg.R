plot.cv.l2.reg <- function(x, ...)
{
 plot.args<-list(log(x$lam.vec),x$mean.error,xlab="log(Lambda)",ylab="CV Error")
  
  do.call("plot", plot.args) 
  do.call("lines", plot.args)
}