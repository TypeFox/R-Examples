plot.quantileDA<-function(x,...)
{
	
	out=x
	if (is.null(out$me.median)) {out.m<-theta.cl(out$train,out$test,out$cl.train.0,0.5,out$cl.test.0)
		mediana<-out.m$me.test} else mediana<-mean(out$me.median)
	
	
	if (is.null(out$me.centroid)) {out.c<-centroidcl(out$train,out$test,out$cl.train.0,out$cl.test.0)
		centroid<-out.c$me.test} else centroid<-mean(out$me.centroid)
	
	max.sup<-max(out$test.rates)+0.4
	plot(out$thetas,out$test.rates,type="l",ylab="misclassification rates",xlab=expression(theta),ylim=c(0,max.sup))
	
	abline(h=centroid,col=2,lty=3)
	abline(h=mediana,col=2,lty=1)
	abline(h=mean(out$me.test),col=3)
	legend(0.7,max.sup,c("centroid","median","best quantile"),col=c(2,2,3),lty=c(3,1,1),bty="n")
}
