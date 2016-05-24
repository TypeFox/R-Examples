# The 'stcor.test' function is an extension of the Box-Pierce statistic test for
# univariate time series, using the results of (Pfeifer & Deutsch, 1981).

stcor.test <- function(data, wlist, tlag=NULL, slag=NULL, fitdf=0) UseMethod("stcor.test")

stcor.test.default <- function(data, wlist, tlag=NULL, slag=NULL, fitdf=0) {

	# Default arguments
	if (is.null(tlag))
		tlag <- floor(10 * log10(nrow(data)))
	if (is.null(slag))
		slag <- length(wlist) - 1

	# Compute the statistic
	cor <- stacf(data, wlist[1:(slag+1)], tlag.max=tlag, plot=F)
	fac <- nrow(data) - 1:tlag %*% t(rep(1, slag+1))
	df <- tlag * (slag+1) - fitdf

	statistic <- ncol(data) * sum( fac * cor^2 )
	pval <- 1 - pchisq(statistic, df)

	out <- data.frame("X-squared"=statistic,
				df=df,
				p.value=pval)

	class(out) <- "stcor.test"
	out

}

print.stcor.test <- function(x, ...) {

	# Print a nice summary
	cat("\tMultivariate Box-Pierce Non Correlation Test\n")
	cat("\t--------------------------------------------\n\n")

	class(x) <- "data.frame"
	print.data.frame(x)

	cat("\nDecision: ")

	if (x$p.value < .05)
		cat("Non Correlation Hypothesis should be rejected.\n")
	else
		cat("Can't reject Non Correlation Hypothesis.\n")

}