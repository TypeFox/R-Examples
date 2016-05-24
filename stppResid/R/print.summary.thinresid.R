print.summary.thinresid <- function(x, ...)
{
	cat("Thinning rate: ", x$k, "\n")
	cat("Number of residuals: ", x$n, "\n")
	cat("Expected number of residuals: ", x$n.exp, "\n")
	cat("One-tailed p-value: ", x$p.val)
}