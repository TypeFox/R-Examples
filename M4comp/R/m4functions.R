#' Plot an M4 time series
#'
#' \code{plot.M4ts} plots a time series from the M4 competition data set.
#' @param x an M4ts object
#' @param xlim Limits on x-axis
#' @param ylim Limits on y-axis
#' @param main Main title
#' @param xlab Label on x-axis
#' @param ylab Label on y-axis
#' @param ...  Other plotting arguments
#' @return None. The function produces a time series plot of the M4ts object \code{x}.
#' @examples
#' plot(M4[[1]])
#' @export

plot.M4ts <- function(x, xlim = c(tsp(x$past)[1], tsp(x$future)[2]), ylim = range(x$past, x$future), main = x$id, xlab = "", ylab = x$units, ...){
    freq <- frequency(x$past)
    plot(ts(c(x$past, rep(NA, x$H)), end = tsp(x$past)[2] + x$H/freq, frequency = freq), ylim = ylim,
                xlim = xlim, ylab = "", xlab = "", ...)
    
    title(main = list(main, cex = 1, font = 2))
    
    lines(ts(x$future, start = tsp(x$past)[2] + 1/freq, frequency = freq), lt = 1, col = 2)
}

#' @export
print.M4ts <- function(x, ...){
	cat(paste("ID    : ", x$id,     "\n"))
	cat(paste("Type  : ", x$type,   "\n"))
	cat(paste("Period: ", x$period, "\n"))
	cat(paste("Units: ", x$units,   "\n\n"))
	
	cat("HISTORICAL: \n")
    print(x$past)
    cat("\nFUTURE:     \n")
    print(x$future)
}

#' @export
print.M4data <- function(x, ...){
	allperiod <- factor(sapply(x, "[[", "period"), levels = c("YEARLY", "BIANNUALLY", "QUARTERLY", "MONTHLY", "WEEKLY", "DAILY"))
	alltypes  <- factor(sapply(x, "[[", "type")  , levels = c("BUSINESS-INDUSTRY", "CLIMATE", "DEMOGRAPHICS", "ECONOMICS", "FINANCE", "INTERNET-TELECOM", "INVENTORY"))
	print(table(alltypes, allperiod, dnn = c("Type", "Period")))
}
