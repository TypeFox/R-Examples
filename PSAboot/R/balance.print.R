#' Print method for balance.
#' 
#' This is a crude measure of overall balance. Absolute value of the standardized
#' effect sizes are calculated for each covariate. Overall balance statistics are
#' the mean of those effect sizes after adjustement for each method across all
#' bootstrap samples.
#' 
#' @param x results from \code{\link{balance}}.
#' @param na.rm whether NA balance statistics should be removed before averaging them.
#' @param ... currently unused.
#' @method print PSAboot.balance
#' @export
print.PSAboot.balance <- function(x, na.rm=TRUE, ...) {
	cat(paste0('Unadjusted balance: ', x$pool.fun(x$unadjusted), '\n'))
	print(cbind(Complete=apply(x$complete, 1, x$pool.fun, na.rm=na.rm),
				Bootstrap=apply(x$pooled, 2, x$pool.fun)))
}
