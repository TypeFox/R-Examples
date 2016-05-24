print.CADFtestsummary <- function(x, ...)
{
	# x is an object of class `CADFtestsummary'
	ttype <- "Covariate Augmented DF test"
	if (nrow(x$test.summary)==3) ttype <- "Augmented DF test"
	cat(ttype, "\n")
	print(x$test.summary, ...)
	print(x$model.summary, ...)
}
