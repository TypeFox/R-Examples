print.structure <-
function(x,...) 
{
	cat("\nLikelihood Ratio Test:\n")
	cat("NH: structure=", colnames(x$AIC)[1], "\n")
	cat("AH: structure=", colnames(x$AIC)[2], "\n")

	print(x$LRT)
	cat("\nInformation Criteria:\n")

	outic <- data.frame(rbind(x$AIC, x$BIC))
  rownames(outic) <- c("aic", "bic")
	print(outic)
	invisible(x)
}
