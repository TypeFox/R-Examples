#' @name chisquare
#' @title Chi-square
#' @param x  input data frame
#' @keywords internal

## Hardcoding for column order should be removed from these.

chisquare <- function(x)
{
	outputmatrix <- matrix(nrow=ncol(x)-1, ncol=2)
	outputmatrix[1,1] <- "Response Metric"
	outputmatrix[1,2] <- "Chi-Square Test P-value"

	iterations <- ncol(x)-2
	for (i in 1:iterations)
	{
	lm <- stats::lm(x[,i+2] ~ x[,2], data=x)
 	outputmatrix[i+1,1] <- colnames(x[i+2])
	bartlett_p <- stats::bartlett.test(x[,i+2] ~ x[,2], data=x)$p.value
	shapiro_p <- stats::shapiro.test(lm$residuals)$p.value
 	outputmatrix[i+1,2] <- stats::pchisq(c((-2*(log(bartlett_p)+log(shapiro_p)))), df=4, lower.tail=FALSE)
	}
	return(outputmatrix)
}
