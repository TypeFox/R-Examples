expectedplot <-
function(OBJ, axis, quants, main, xlim, ylim, xlab, ylab, ...){


		if(missing(main)){main="Expected Total Score\n"}
		if(missing(ylab)){ylab="Expected Score"}
		if(missing(ylim)){ylim=c(min(OBJ$subjscore),max(OBJ$subjscore))}
		if(missing(xlim)){xlim=c(min(axis),max(axis))}	

	
	plot(axis,OBJ$expectedscores,type="l",ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim,main=main,...)

	
	axis(3,at=quants, lab=labels(quants),tck=0)
	abline(v=quants,col="blue",lty=2)
	
	return(round(OBJ$expectedscores,3))
}

