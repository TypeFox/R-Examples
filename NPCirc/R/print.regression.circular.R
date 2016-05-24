print.regression.circular<-function (x, digits = NULL, ...) {
	cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
	" (", x$n, " obs.);", "\tBandwidth 'bw' = ", formatC(x$bw, 
	digits = digits), "\n\n", sep = "")
	print(summary(as.data.frame(x[c("x", "y")])), digits = digits, ...)
	invisible(x)
}
