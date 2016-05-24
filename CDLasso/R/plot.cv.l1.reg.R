plot.cv.l1.reg <- function(x, ...)
{
 plot.args<-list(x$lam.vec,x$mean.error,xlab="Lambda",ylab="CV Error")
  
  do.call("plot", plot.args) 
  do.call("lines", plot.args)
}