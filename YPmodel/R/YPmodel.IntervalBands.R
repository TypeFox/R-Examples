YPmodel.IntervalBands <-
function(data=c(), Internal=c(), Estimate=c(), ...)
{


	if(is.null(data)){
		stop(paste(fun.errorMessage('DataSet')))
	}

	Data <- YPmodel.inputData(data=data)

#	if(is.null(Parameters)){
#			warning(paste(fun.errorMessage('DefaultParameter')))
#			Parameters <- YPmodel.setParameter(Data=Data)
#	}

	if(is.null(Estimate)){
		Estimate <- YPmodel.estimate(data=data)
	}


	if(is.null(Internal)){
		Internal <- fun.internalParameters(Data=Data, Estimate=Estimate)
	}
	h <- 0.001

	ru <- Internal$ru
	p <- Internal$p
	pl <- Internal$pl
	bt <- Internal$bt
	deni <- Internal$deni
	kall <- Internal$kall
	sm <- Internal$sm
	gama <- Internal$gama
	b <- Internal$b
	GroupData <- Internal$GroupData

	data8 <- fun.hazardRationEstimation(Data,GroupData,ru,p,pl,bt,deni,kall,sm,gama,b,h)
	hr <- data8$hr
	ld2 <- data8$ld2
	ud2 <- data8$ud2
	upp22 <- data8$upp22
	low22 <- data8$low22
	upp3 <- data8$upp3
	low3 <- data8$low3
	upp90 <- data8$upp90
	low90 <- data8$low90

	Result <- c()
	Result$hr <- hr
	Result$ld2 <- ld2
	Result$ud2 <- ud2
	Result$upp22 <- upp22
	Result$low22 <- low22
	Result$upp3 <- upp3
	Result$low3 <- low3
	Result$upp90 <- upp90
	Result$low90 <- low90

	IntervalBands <- Result
	IntervalBands$Data <-Data
	IntervalBands$Estimate <- Estimate

	class(IntervalBands) <- "YPmodel.IntervalBands"
	IntervalBands$call <- match.call()
	#-----------------------------------------------------------------#
	## Output Resuts
	#-----------------------------------------------------------------#
	    return(IntervalBands)
	#-----------------------------------------------------------------#


}
