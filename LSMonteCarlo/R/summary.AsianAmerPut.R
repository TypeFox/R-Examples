summary.AsianAmerPut <-
function(object, ...){
	cat("\nAsian American Put Option\n(Pay-off with arithmetic mean)\nMethod: Simple Least Squares Monte Carlo\n")
	mat1<-matrix(object)
	rownames(mat1)<-c("Option price", "Spot price", "Strike", "Volatility", "Number of paths", "Number of time-steps", "Interest rate", "Dividend rate", "Maturity time")
	colnames(mat1) <- c(" ")
	print(mat1)
}
