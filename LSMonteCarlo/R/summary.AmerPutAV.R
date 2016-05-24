summary.AmerPutAV <-
function(object, ...){
	cat("\nAmerican Put Option\nMethod: Least Squares Monte Carlo with Antithetic Variates\n")
	mat1<-matrix(object)
	rownames(mat1)<-c("Option price", "Spot price", "Strike", "Volatility", "Number of original paths", "Number of antithetic paths", "Number of time-steps", "Interest rate", "Dividend rate", "Maturity time")
	colnames(mat1) <- c(" ")
	print(mat1)
}
