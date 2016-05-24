summary.YPmodel.estimate <-
function(object=c(), ...)
{
	Data <- object$Data
	interval <- object$interval

	Estimate <- object

	Beta <- Estimate$beta
	beta <- Estimate$beta
	n <- Data$length
	variance.beta1 <- Estimate$variance.beta1
	variance.beta2 <- Estimate$variance.beta2

	colnames(Beta) <- c("Beta_1","Beta_2")
	rownames(Beta) <- c("estimates")
	BetaInterval <- matrix(c(beta[1,1]-1.96*(variance.beta1/sqrt(n)),
		beta[1,2]-1.96*(variance.beta2/sqrt(n)),
		beta[1,1]+1.96*(variance.beta1/sqrt(n)),
		beta[1,2]+1.96*(variance.beta2/sqrt(n))
		), nrow=2, ncol=2)
	colnames(BetaInterval) <- c("Lower bound", "Upper bound")
	rownames(BetaInterval) <- c("Beta_1","Beta_2")

	Theta <- matrix(exp(beta), nrow=1, ncol=2)
	colnames(Theta) <- c("Theta_1","Theta_2")
	rownames(Theta) <- c("estimates")
	ThetaInterval<- exp(BetaInterval)
	colnames(ThetaInterval) <- c("Lower bound", "Upper bound")
	rownames(ThetaInterval) <- c("Theta_1","Theta_2")


	cat("\n-------------------------------------------------------------------------------------------------------------  \n")
	cat("Parameters of short-term and long-term hazard ration model  \n")
	cat("\n Adaptive weight (Beta):\n \n")
	printCoefmat(Beta, digits=4)
	cat("\n Hazard ratio:\n \n")
	printCoefmat(Theta)
	if(interval==1){
		cat("\n Confidence Interval (Beta):\n \n")
		printCoefmat(BetaInterval)
		cat("\n Confidence Interval (Theta)\n \n")
		printCoefmat(ThetaInterval)
	}

}
