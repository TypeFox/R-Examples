summary.AmerPutCV <-
function(object, ...){
	cat("\nAmerican Put Option\nMethod: Least Squares Monte Carlo with Control Variates\n(B&S solution for European Put is used as a control)\n")
	mat1<-matrix(object)
	rownames(mat1)<-c("Option price", "Spot price", "Strike", "Volatility", "Number of paths", "Number of time-steps", "Interest rate", "Dividend rate", "Maturity time")
	colnames(mat1) <- c(" ")
	print(mat1)
}
