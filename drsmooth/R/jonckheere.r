#' @name jonckheere
#' @title Jonckheere
#' @param x input data frame
#' @keywords internal

## Hardcoding for column order should be removed from these.
## Also, this code appears to assume that dose-response pairs are already in rank order
## Unclear how ties would be handled.

jonckheere <- function(x)
{
	outputmatrix <- matrix(nrow=ncol(x)-1, ncol=2)
	outputmatrix[1,1] <- "Response Metric"
	outputmatrix[1,2] <- "Jonckheere Test Trend & P-value"

	iterations <- ncol(x)-2
	for (i in 1:iterations)
	{
	lm <- lm(x[,i+2] ~ x[,2], data=x)
 	outputmatrix[i+1,1] <- colnames(x[i+2])
  Cor_test <- stats::cor.test(x[,i+2], x[,1], method="k")
  if (Cor_test$statistic > 0) {Trend <- "Positive"} else {Trend <- "Negative"}
 	outputmatrix[i+1,2] <- paste(Trend, Cor_test$p.value, sep = "         ")
	}
	return(outputmatrix)
}