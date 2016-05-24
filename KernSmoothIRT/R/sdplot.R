sdplot <-
function(OBJ, quants, axis, main, xlab,ylab,...){

	
	Testvar <- apply(OBJ$OCC[,-c(1:3)],2,function(x)sum(x*OBJ$OCC[,3]**2 - (x*OBJ$OCC[,3])**2))



	Testsd<-sqrt(Testvar)

	if(missing(main)){main="Test Standard Deviation\n"}
	if(missing(ylab)){ylab="Standard Deviations"}

	plot(axis,Testsd,type="l",main=main,xlab=xlab,ylab=ylab,...)

		axis(3,at=quants, lab=labels(quants),tck=0)
		abline(v=quants,col="blue",lty=2)
	

	return(round(Testsd,3))	
}

