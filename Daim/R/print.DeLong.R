
print.DeLong <- function(x, digits=max(3, getOption("digits") - 3), ...){
	cat("Estimated AUC's:\n")
	print(format(round(x$AUC, digits=digits, ...), nsmall=digits, ...))
	cat("\n Pairwise comparisons:\n")
	print(format(round(x$difference, digits=digits, ...), nsmall=digits, ...))
	cat(paste("\n Overall test:\n p-value =", format.pval(x$global.p, digits=digits), "\n"))
	invisible(x)
}

