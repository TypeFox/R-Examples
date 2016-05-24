print.sglasso <- function (x, digits = max(3, getOption("digits") - 3), ...){
	rho <- x$rho
	df <- x$df
	x$algorithm
	tbl <- data.frame(rho, df)
	names(tbl) <- c("rho", "df")
	tbl.format <- format(tbl, digits = digits)
	cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	print(tbl.format, print.gap = 2, quote = FALSE, row.names = FALSE, ...)
	cat("\nAlgorithm", x$algorithm, "with exit =", x$conv, "\n\n")
	invisible(tbl)
}