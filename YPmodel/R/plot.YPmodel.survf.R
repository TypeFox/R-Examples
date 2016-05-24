plot.YPmodel.survf <-
function(x, Internal=c(), ...)
{



	Data <- x$Data
	rTestData <- x$rTestData
	Estimate <- x$Estimate

	LackFitTest <- x

	if(is.null(Data)){
		stop(paste(fun.errorMessage('DataSet')))
	}

	if(is.null(Internal)){
		Internal <- fun.internalParameters(Data=Data, Estimate=Estimate)
	}


	X <- Data$X
	kk <- Internal$kk
	wtildCount2 <- LackFitTest$wtildCount2
	lineCount2 <- LackFitTest$lineCount2
	obs2 <- LackFitTest$obs2

	#dev.new()
	plot(X[1:kk]*365,wtildCount2[,1],"l",lty = "dotted",col="red",xlab="Days", ylab=" ",xlim=c(1, max(X*365)),ylim=c(min(obs2,wtildCount2), max(obs2,wtildCount2)))
	for(i in 2:lineCount2){
	    lines(X[1:kk]*365,wtildCount2[,i],"l",lty = "dotted",col="red")
	}
	lines(X[1:kk]*365,obs2[1:kk],"l",col="blue",lwd=2)
	title(main="Plots of the contrast-based test statistic \r\n and randomly selected realizations of the process")


}
