print.inflmlm <-
function(x, digits = max(3, getOption("digits") - 4), FUN=det, ...) {
	df <- as.data.frame(x, FUN=FUN)
	cat("Multivariate influence statistics for model:\n", 
	 paste(deparse(x$call), sep = "\n", collapse = "\n"), 
	 "\n m=", x$m, "case deletion diagnostics",
	 ifelse(x$m>1, paste(", using", deparse(substitute(FUN)), "for matrix values\n\n"), "\n"))
	print(df, digits=digits, ...)
	invisible(x)
}
