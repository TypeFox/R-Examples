#' Lag statistics on timeseries of log returns
#' 
#' Lag statistics on log returns are added to the xts attributes. It takes
#' a vector of lags and calculates the mean, stdev, var, skewness, and kurtosis
#' for cumulative log returns of each lag. The data is stored as a list of vectors 
#' under \code{lagstats} attribute. Be aware this function uses multicore lapply.
#'
#' @param ts the xts object from sample data. The ts must have the logr column.
#'           If a string is given, it will be replaced with sample data of the symbol.
#' @param lags a numeric vector of integers greater than 0.
#' @param absolute logical, if \code{TRUE}, statistics calculated on absolute log returns. 
#'              Default: \code{FALSE}.
#'
#' @return The xts object containing \code{lagstats} attribute
#'
#' @keywords sample-data statistics
#'
#' @export
#'
#' @import zoo xts
#' @importFrom stats var
#'
###  Note: example failed in R CMD check due to missing na.fill()
###' @examples
###' \dontrun{
###' dji <- ecd.ts_lag_stats(ecd.data("dji"), 2)
###' }
### <======================================================================>
ecd.ts_lag_stats <- function(ts = "dji", lags, absolute=FALSE) {

	if (is.character(ts)) {
		ts <- ecd.data_stats(ts) # ts is a symbol, replace it with sample data
	}
    
	logr1d <- ts$logr
	
	calc_lag_stats <- function (lag) {
	    logr <- zoo::rollapply(logr1d, lag, sum, fill=NA)
	    logr <- logr[!is.na(logr$logr)]
        if (absolute) logr <- abs(logr)
        
	    list(
	        mean = unname(as.numeric(mean(logr))),
            var  = unname(as.numeric(var(logr))),
	        skewness = unname(as.numeric(moments::skewness(logr))),
	        kurtosis = unname(as.numeric(moments::kurtosis(logr)))
	    )
	}

	rs <- parallel::mclapply(lags, calc_lag_stats)
    get <- function(s) sapply(rs, function(l) l[[s]])

    # statistics
    V <- get("var")
    lagstats <- list(
        lags     = lags,
        absolute = absolute,
        mean     = get("mean"),
        stdev    = sqrt(V),
        var      = V,
        skewness = get("skewness"),
        kurtosis = get("kurtosis")
    )

    xtsAttributes(ts) <- list(lagstats = lagstats)
    ts
}

