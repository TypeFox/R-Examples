YPmodel.default <-
function(data, startPoint=c(0,0), nm=log(100), maxIter1=50, maxIter2=20, repNum=1000, ...)
{

	Estimate <- YPmodel.estimate(data=data, startPoint=startPoint, nm=nm, maxIter1=maxIter1, maxIter2=maxIter2)

	IntervalBands <- YPmodel.IntervalBands(data=data, Estimate=Estimate) 

	LackFitTest <- YPmodel.lackfittest(data=data, repNum=repNum)

	Adlgrk <- YPmodel.adlgrk(data=data, Estimate=Estimate)

	Result <- c()
	Result$call <- match.call()
	class(Result) <- "YPmodel"
	Result$Data <- YPmodel.inputData(data=data)
	Result$Estimate <- Estimate
	Result$IntervalBands <- IntervalBands
	Result$LackFitTest <- LackFitTest
	Result$Adlgrk <- Adlgrk

	return(Result)
	
}
