summary.QuantoAmerPut_AV <-
function(object, ...){
	cat("\nQuanto American Put Option\nMethod: Least Squares Monte Carlo with Antithetic Variates\n")
	mat1<-matrix(object)
	rownames(mat1)<-c("Option price", "Spot price", "Strike", "Volatility", "Number of original paths", "Number of antithetic paths", "Number of time-steps", "Interest rate", "Dividend rate", "Maturity time","Spot price of 3rd asset", "Volatility of 3rd asset", "Interest rate of 3rd asset", "Dividend rate of 3rd asset","Correlation coefficient")
	colnames(mat1) <- c(" ")
	print(mat1)
}
