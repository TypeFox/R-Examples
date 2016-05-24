error.bars <-
function(x, upper, lower, width = 0.02, ...)
{
	xlim <- range(x)
	barw <- diff(xlim) * width
	segments(x, upper, x, lower, ...)
	segments(x - barw, upper, x + barw, upper, ...)
	segments(x - barw, lower, x + barw, lower, ...)
	range(upper, lower)
}

plot.cv.glmgraph <- function(x,...){

	lambda2 <- x$obj$lambda2
	nlambda2 <- length(x$cvm)
    n <- ceiling(sqrt(nlambda2))    
    lambda1 <- x$obj$lambda1
    type.measure <- x$type.measure
   
   	mypar=function (a = 2, b = 2, brewer.n = 8, brewer.name = "Dark2", ...) {
  		par(mar = c(2.5, 2.5, 1.6, 1.1), mgp = c(1.5, 0.5, 0))
 		par(mfrow = c(a, b), ...)
	}
	mypar(n,n)
    
    for(i in 1:nlambda2){
    	cvm <- x$cvm[[i]]
    	cvsd <- x$cvsd[[i]]
    	lambda1 <- x$obj$lambda1[1:length(cvm)]
    	plot(x=log(lambda1),y=cvm,ylim=range(cvm-cvsd,cvm+cvsd),ylab=type.measure,xaxt = "n",xlab=expression(lambda[1]),type="n",main= substitute(lambda[2]==p, list(p=format(lambda2[i],scientific=TRUE))))
		error.bars(log(lambda1),cvm-cvsd,cvm+cvsd,width=0.01,col="darkgrey")
  		points(log(lambda1),cvm,pch=20,col="red")
		axis(side=1,at=log(lambda1),labels=round(lambda1,digits=3),tick=FALSE,line=0)
		if(type.measure=="auc") lmin=getmin(lambda1,-cvm,cvsd)
		else lmin=getmin(lambda1,cvm,cvsd)
		abline(v=log(lmin$lambda.min),lty=2)
		abline(v=log(lmin$lambda.1se),lty=3)
  		invisible()
    }   
}






















  




