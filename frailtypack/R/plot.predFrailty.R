
"plot.predFrailty" <- function (x, conf.bands=FALSE, pos.legend="topright", cex.legend=0.7, ylim=c(0,1), ...)
{
	if ((conf.bands) & (!x$icproba)) stop("Confidence intervals were not calculated. Use the MC.sample argument in the 'prediction' function")
	
	if (x$moving.window){ 
		legende <- paste("Predicted cumulative probability of event between",x$t,"and time t")
		xlim <- c(x$t-(max(x$x.time)-min(x$x.time))*0.1,max(x$x.time))
	}else{ 
		legende <- paste("Predicted probability of event in the next",x$window)
		xlim <- c(min(x$x.time),max(x$x.time))
	}
	
	# par(mfrow=c(2,1))
	if (is.null(x$type)){
		title <- paste("Dynamic prediction for a Cox model")
	}else{
		if (x$type=="marginal") title <- paste("Marginal dynamic prediction for a shared frailty model")
		else title <- paste("Conditional dynamic prediction for a shared frailty model")
	}
	
	if (conf.bands){
		matplot(x$x.time,t(x$pred),type="l",lty=1,xlab="Time t",ylab=legende,main=title,ylim=ylim,xlim=xlim)
		matlines(x$x.time,t(x$predLow),type="l",lty=2)
		matlines(x$x.time,t(x$predHigh),type="l",lty=2)
	}else{
		matplot(x$x.time,t(x$pred),type="l",lty=1,xlab="Time t",ylab=legende,main=title,ylim=ylim,xlim=xlim)
	}
	
	legend(pos.legend, paste("profile",(1:x$npred)),lty=1,bty="n",col=(1:x$npred))

	if (x$moving.window) abline(v=x$t,lty=2)

	return(invisible())
}
