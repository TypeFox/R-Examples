#' @name bartlett
#' @title Bartlett's
#' @param x  input data frame
#' @keywords internal

## Hardcoding for column order should be removed from these.

bartlett <- function(x)
{
	outputmatrix <- matrix(nrow=ncol(x)-1, ncol=2)
	outputmatrix[1,1] <- "Response Metric"
	outputmatrix[1,2] <- "Bartlett Test P-value"

	iterations <- ncol(x)-2
	for (i in 1:iterations)
	{
	lm <- lm(x[,i+2] ~ x[,2], data=x)
 	outputmatrix[i+1,1] <- colnames(x[i+2])
 	outputmatrix[i+1,2] <- bartlett_p <- stats::bartlett.test(x[,i+2] ~ x[,2], data=x)$p.value
	}
	return(outputmatrix)
}
