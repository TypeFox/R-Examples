

plot.Daim.vector <- function(x, color="blue", type="l", bty="n", 
		xlab="False positive rate",ylab="True positive rate", main="ROC curve", ...)
{
	best.cut <- which.max((1-x$FPR)+x$TPR-1)
	plot(x$FPR, x$TPR, xlab=xlab, ylab=ylab, main=main, col=color, type=type, ...)
	points(x$FPR[best.cut], x$TPR[best.cut], col=1, pch=19)
	grid()
    lines(c(0,1),c(0,1),col="black")
	legend("bottomright",legend=
		c(paste("best cut-off:",formatC(x$cutoff[best.cut],digits=max(3, getOption("digits") - 3))),
		paste("FPR:",formatC(x$FPR[best.cut],digits=max(3, getOption("digits") - 3))),
		paste("TPR:",formatC(x$TPR[best.cut],digits=max(3, getOption("digits") - 3))),
		paste("AUC:",formatC(auc(1-x$FPR,x$TPR),digits=max(3, getOption("digits") - 3)))),inset=0.01,...)
}



plot.Daim.list <- function(x, color=rgb(1,0,0,alpha=0.5),
			lty=1, lwd=1, pch=19,
			xlab="False positive rate",ylab="True positive rate", 
			main="ROC curves", legend=TRUE, ...)
{
	plot(c(0,1),c(0,1), type="n",xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, main=main)
	grid()
    lines(c(0,1),c(0,1),col="black")
    tmp.col <- color
	tmp.lty <- lty
	tmp.lwd <- lwd
	tmp.pch <- pch
	for(i in 1:length(x)){
		if(length(color) > 1 && length(color) >= i)
		tmp.col <- color[i]
		if(length(lty) > 1 && length(lty) >= i)
		tmp.lty <- lty[i]
		if(length(lwd) > 1 && length(lwd) >= i)
		tmp.lwd <- lwd[i]
		if(length(pch) > 1 && length(pch) >= i)
		tmp.pch <- pch[i]
		best.cut <- which.max((1-x[[i]]$FPR)+x[[i]]$TPR-1)
		lines(x[[i]]$FPR, x[[i]]$TPR, col=tmp.col, lty=tmp.lty, lwd=tmp.lwd, ...)
		points(x[[i]]$FPR[best.cut], x[[i]]$TPR[best.cut], col=tmp.col, pch=tmp.pch, ...)
	}
	if(legend && length(x) < 11){
		Nvar <- names(x)
		legend("bottomright", paste("best cut-off ",Nvar,sep=""),  col=color,  
			   lty=lty,pch=pch, merge=TRUE, inset=0.01, bg="white", ...)
	}
}