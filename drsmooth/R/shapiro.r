#' @name shapiro
#' @title Shapiro
#' @param x  input data frame
#' @keywords internal

## Hardcoding for column order should be removed from these.

shapiro <- function(x)
{
	outputmatrix <- matrix(nrow=ncol(x)-1, ncol=2)
	outputmatrix[1,1] <- "Response Metric"
	outputmatrix[1,2] <- "Shapiro Test P-value"

	iterations <- ncol(x)-2
	for (i in 1:iterations)
	{
	lm <- stats::lm(x[,i+2] ~ x[,2], data=x)
 	outputmatrix[i+1,1] <- colnames(x[i+2])
 	outputmatrix[i+1,2] <- shapiro_p <- stats::shapiro.test(lm$residuals)$p.value
	}
	return(outputmatrix)
}
