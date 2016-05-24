plot.YPmodel.martint <-
function(x=c(), Internal=c(), ...)
{

	Data <- x$Data
	rTestData <- x$rTestData
	Estimate <- x$Estimate

	LackFitTest <- x

	#-----------------------------------------------------------------#

	if(is.null(Internal)){
		Internal <- fun.internalParameters(Data=Data, Estimate=Estimate)
	}

	X <- Data$X
	kk <- Internal$kk
	wtildCount1 <- LackFitTest$wtildCount1
	lineCount1 <- LackFitTest$lineCount1
	obs <- LackFitTest$obs

	#dev.new()
	plot(X[1:kk]*365,wtildCount1[,1],"l",lty = "dotted",col="red",xlab="Days", ylab=" ",xlim=c(1, max(X*365)),ylim=c(min(obs,wtildCount1), max(obs,wtildCount1)))
	for(i in 2:lineCount1){
	    lines(X[1:kk]*365,wtildCount1[,i],"l",lty = "dotted",col="red")
	}
	lines(X[1:kk]*365,obs[1:kk],"l",col="blue",lwd=2)
	title(main="Plots of the martingale residual-based test statistic \r\n and randomly selected realizations of the process")
}
