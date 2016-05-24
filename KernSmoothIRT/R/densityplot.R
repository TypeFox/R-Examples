densityplot <-
function(x,xlim,ylim,xlab,ylab,main,...){



		if(missing(main)){main="Observed Score Distribution\n"}
	#	if(missing(xlab)){xlab="Score"}
		if(missing(xlim)){xlim=c(0,max(x$subjscore))}
		
		if(missing(ylab)){ylab="Density of Score"}
	

		xlab <- "Observed Scores"

		ymax<-max(density(x$subjscore,from=0,to=max(x$subjscore))$y)


		if(missing(ylim)){ylim=c(0,ymax)}
	
	plot(density(x$subjscore,from=0,to=max(x$subjscore)),xlim=xlim,ylim=ylim,ylab=ylab,xlab=xlab,main = main,...)

		axis(3,at=x$subjscoresummary, labels =labels(x$subjscoresummary),tck=0)
		abline(v=x$subjscoresummary,col="blue",lty=2)
		box()

}

