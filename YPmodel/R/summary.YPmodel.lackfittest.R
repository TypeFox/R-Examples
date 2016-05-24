summary.YPmodel.lackfittest <-
function(object=c(),...)
{

	LackFitTest <- object

	newBeta <- LackFitTest$newBest
	pvalu1 <- LackFitTest$pvalu1
	pvalu2 <- LackFitTest$pvalu2

	colnames(newBeta) <- c("Beta_1","Beta_2")
	rownames(newBeta) <- c("estimates")


	cat("\n-------------------------------------------------------------------------------------------------------------  \n")
	cat("Lack-of-fit tests for checking short-term and long-term hazard ration model  \n")
	cat("\n Adaptive weight (Beta, sample odds function estimator using only the control group data): \n \n")
	printCoefmat(newBeta, digits=4)
	cat("\n Residual, the martingale residual-based test (p-value):\n")
	cat(paste(pvalu1),"\n")
	cat("\n Contrast, the contrast-based test (p-value):\n")
	cat(paste(pvalu2),"\n")

}
