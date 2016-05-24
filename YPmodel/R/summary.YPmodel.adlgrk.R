summary.YPmodel.adlgrk <-
function(object=c(), ...)
{
	Adlgrk <- object

	pval <- Adlgrk$pval

	cat("-------------------------------------------------------------------------------------------------------------  \n")
	cat("\nImproved Logrank-Type Tests (p-value):\n")
	cat(paste(round(pval,digits = 4)),"\n")

}
