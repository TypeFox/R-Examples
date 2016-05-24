plot.survAUC <- function(x, col="red", type="l", ylim=c(0,1),
					xlab="Time",ylab="AUC", main="Time-dependent AUC", add=FALSE, ...)
{
	if(!add){
		plot(x$times, x$auc, ylim=ylim, main=main, xlab=xlab, ylab=ylab, col=col, type=type, ...)
	}
	else{
		lines(x$times, x$auc, col=col, ...)
	}
}





plot.survErr <- function(x, col="red", type="l", ylim=c(0,1),
xlab="Time",ylab="Prediction error", main="Time-dependent Prediction Error", add=FALSE, ...)
{
	if(!add){
		plot(x$times, x$error, ylim=ylim, main=main, xlab=xlab, ylab=ylab, col=col, type=type, ...)
	}
	else{
		lines(x$times, x$error, col=col, ...)
	}
}

