plot.cv.vda.r <- function(x, ...)
{
 plot.args<-list(x$lam.vec,x$mean.error,xlab="log(Lambda)",ylab="CV Error")
  
  do.call("plot", plot.args) 
  do.call("lines", plot.args)
}
