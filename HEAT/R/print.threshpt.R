print.threshpt <-
function(x, ...)
{
	cat("Call:\n")
	print(x$call)
	cat("\nCoefficietns for main exposure:\n")
	print(x$best.fit[1:6])
	cat("\nOther coefficients:\n")
	print(x$parm.coef)
}
