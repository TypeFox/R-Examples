YPmodel.lackfittest <-
function(data=c(), repNum=1000, Internal=c(), ...)
{
	if(is.null(data)){
		stop(paste(fun.errorMessage('DataSet')))
	}

	Data <- YPmodel.inputData(data=data)


	rTestData <- YPmodel.setRandom(Data=Data, repNum=repNum)


	if(is.null(Internal)){
		Estimate <- YPmodel.estimate(data=data)
		Internal <- fun.internalParameters(Data=Data, Estimate=Estimate)
	}

	jh <- 0.001

	u0 <- Internal$u0
	fb <- Internal$fb
	kall <- Internal$kall
	kk <- Internal$kk
	l <- Internal$l
	dfb1 <- Internal$dfb1

	data11 <- fun.newBest(Data,u0)
	newBest <-data11$best

	data9 <- fun.martint(Data,newBest,u0,fb,kall,kk,rTestData, repNum, jh)
	mobs1 <- data9$mobs1
	pvalu1 <- data9$pvalu1
	obs <- data9$obs
	wtildCount1 <- data9$wtildCount
	lineCount1 <- data9$lineCount

	data10 <- fun.survf(Data,newBest,u0,fb,kall,kk,l,dfb1,rTestData, repNum, jh)
	mobs2 <- data10$mobs2
	pvalu2 <- data10$pvalu2
	obs2 <- data10$obs2
	wtildCount2 <- data10$wtildCount
	lineCount2 <- data10$lineCount

	Result <- c()
	Result$newBest <- newBest
	Result$pvalu1 <- pvalu1
	Result$pvalu2 <- pvalu2
	Result$mobs1 <- mobs1
	Result$mobs2 <- mobs2
	Result$obs <- obs
	Result$obs2 <- obs2

	Result$wtildCount1 <- wtildCount1
	Result$lineCount1 <- lineCount1
	Result$wtildCount2 <- wtildCount2
	Result$lineCount2 <- lineCount2

	LackFitTest <- Result
	LackFitTest$Data <- Data
	LackFitTest$rTestData <- rTestData
	LackFitTest$Estimate <- Estimate

	class(LackFitTest) <- "YPmodel.lackfittest"
	LackFitTest$call <- match.call()

	return(LackFitTest)

}
