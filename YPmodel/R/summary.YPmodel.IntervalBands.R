summary.YPmodel.IntervalBands <-
function(object=c(), Internal=c(), ...)
{
	Data <- object$Data 
	Estimate <- object$Estimate 
	IntervalBands <- object

	if(is.null(Internal)){
		Internal <- fun.internalParameters(Data=Data, Estimate=Estimate)
	}

	X <- Data$X
	ld2 <- IntervalBands$ld2
	ud2 <- IntervalBands$ud2
	hr <- IntervalBands$hr
	upp3 <- IntervalBands$upp3
	low3 <- IntervalBands$low3
	upp22 <- IntervalBands$upp22
	low22 <- IntervalBands$low22
	upp90 <- IntervalBands$upp90
	low90 <- IntervalBands$low90


	Days <- X[ld2:ud2]*365
	hazardRationFunction <- hr[ld2:ud2]
	PointwiseConfidenceLimitsU <- upp3[ld2:ud2]
	PointwiseConfidenceLimitsL <- low3[ld2:ud2]
	ConfidenceBandsU95 <- upp22[ld2:ud2]
	ConfidenceBandsL95 <- low22[ld2:ud2]
	ConfidenceBandsU90 <- upp90[ld2:ud2]
	ConfidenceBandsL90 <- low90[ld2:ud2]
	IntervalsPrint <- cbind(Days, hazardRationFunction, PointwiseConfidenceLimitsU,PointwiseConfidenceLimitsL,ConfidenceBandsU95,ConfidenceBandsL95,ConfidenceBandsU90,ConfidenceBandsL90)
	colnames(IntervalsPrint) <- c("Days", "HR_fun", "lower.cl", "upper.cl", "lower.95%band", "upper.95%band","lower.90%band", "upper.90%band")

	cat("\n-------------------------------------------------------------------------------------------------------------  \n")
	cat("Point estimates, Pointwise confidence intervals, and confidence bands of short-term and long-term hazard ration model  \n\n")
	printCoefmat(IntervalsPrint, digits=4)

}
