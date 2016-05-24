
"plot.epoce" <- function (x, type, pos.legend="topright", cex.legend=0.7, ...)
{

	if (x$new.data){
    if(!missing(type) && type!="mpol")stop("For epoce with new dataset only mpol is calculated and can be plotted")
		plot(x$pred.times,x$mpol,type="b",pch="X",col="blue",xlab="time",ylab="epoce")
		legend(pos.legend,c("mpol"),lty=1,col="blue")
	}else{
    if(!missing(type) && !type%in%c("cvpol","mpol"))stop("Only mpol or cvpol can be plotted for the epoce estimator")
		plot(x$pred.times,noquote(paste("x$",type,sep="")),type="b",pch="X",col="red",xlab="time",ylab="epoce")
		legend(pos.legend,type,lty=1,col="red")
	}

	return(invisible())
}
