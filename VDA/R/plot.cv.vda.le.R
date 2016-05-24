plot.cv.vda.le <-
function(x, ...)
{
 if (sum(x$lam.vec.1)==0){
 	 plot.args<-list(log(x$lam.vec.2),x$error.cv, xlab="log(Lambda2)", ylab="CV Error", col="red")
  	do.call("plot", plot.args) 
  	do.call("lines", plot.args)
 }
 else if (sum(x$lam.vec.2)==0){
 	 plot.args<-list(log(x$lam.vec.1),x$error.cv, xlab="log(Lambda1)", ylab="CV Error", col="red")
	do.call("plot", plot.args) 
  	do.call("lines", plot.args)
 }
 else{
 plot.args1<-list(log(x$lam.vec.1),log(x$lam.vec.2),x$error.cv, xlab="log(Lambda1)", ylab="log(Lambda2)", zlab="CV Error", col="red")
 
 plot.args2<-list(log(x$lam.vec.1),log(x$lam.vec.2),x$error.cv, xlab="log(Lambda1)", ylab="log(Lambda2)", zlab="CV Error", col="blue")
   
  do.call("plot3d", plot.args1) 
  do.call("surface3d", plot.args2)
 } 
}
