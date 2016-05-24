YPmodel.setRandom <-
function(n=c(), repNum=c(), Data=c(), ...)
{

	if(is.null(n)){
		if(! is.null(Data)){
			n <- Data$length
		}
		else{
			stop(paste(fun.errorMessage('DataLength')))
		}
	}

	rTestData <- c()
	rTestData$gaussianRandomDataAll <- rnorm(n*repNum, mean = 0, sd = 1)
	rTestData$uniformRandomData1 <- runif(2, min = 0, max = 1) 
	rTestData$uniformRandomData2 <- runif(100, min = 0, max = 1) 
	rTestData$repNum <- repNum

	return(rTestData)

}
