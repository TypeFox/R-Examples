summary.PriceSurface <-
function(object, ...){
	cat("\nVolatility sequence:\n")
	cat(as.numeric(rownames(object)))
	cat("\n\nStrike sequence:\n")
	cat(as.numeric(colnames(object)))
	cat("\n\nAverage price: ")
	cat(mean(object))
	cat("\nMinimum price: ")
	cat(min(object))
	cat("\nMaximum price: ")
	cat(max(object))
}
