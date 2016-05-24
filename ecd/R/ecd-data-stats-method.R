#' Statistics and histogram on log returns
#' 
#' Statistics and histogram on log returns are added to the xts attributes
#'
#' @param ts can be either a symbol of sample data,
#'           or the xts object from sample data
#' @param breaks A length-one numeric, breaks for generating the histogram.
#' @param merge_tails A length-two numeric vector.
#'        The first element is how many points in the left tail of histogram to be dropped during fitting.
#'        The second element is how many points in the right tail of histogram to be dropped during fitting.
#' @param with.tail logical, include tail statistics, mainly on asypmtotic kurtosis. Default: \code{FALSE}.
#' @param tail.N1 a numeric, defining the wider range of tail statistics
#' @param tail.N2 a numeric, defining the smaller range of tail statistics
#'
#' @return The xts object containing ecd added attributes
#'
#' @keywords sample-data statistics
#'
#' @export
#'
#' @importFrom utils head
#' @importFrom stats var lm predict
#'
#'
#' @examples
#' dji <- ecd.data_stats(ecd.data("dji"))
#' dji <- ecd.data_stats("dji")
### <======================================================================>
ecd.data_stats <- function(ts = "dji", breaks = 20, merge_tails = c(0,0),
                           with.tail=FALSE, tail.N1=7, tail.N2=5) {

	if (is.character(ts)) {
		ts <- ecd.data(ts) # ts is a symbol, replace it with sample data
	}
    
    logr <- as.vector(ts$logr)
	tail_quantile <- 1/length(logr)
        
    # statistics
    v <- unname(as.numeric(var(logr)))
    xtsAttributes(ts) <- list(
        mean     = unname(as.numeric(mean(logr))),
        stdev    = sqrt(v),
        var      = v,
        skewness = unname(as.numeric(moments::skewness(logr))),
        kurtosis = unname(as.numeric(moments::kurtosis(logr))),
        tail_quantile = tail_quantile
    )

    # histogram
    h <- hist(ts$logr, breaks = breaks, plot = FALSE)
    xtsAttributes(ts) <- list(hist = h, breaks = breaks)

	# histuple is an internal representation of histogram
    htu0 <- list(hx = h$mids, hy = h$counts)
    
    htu <- ecd.manage_hist_tails(htu0, merge_tails)
    xtsAttributes(ts) <- list(histuple = htu, merge_tails = merge_tails)

    if (with.tail) {
        n <- seq(1,2^tail.N1)
        abs_logr <- rev(sort(abs(logr)))
        stats_n <- function(n) {
            cutoff <- 999.0
            if ( n > 0 ) {
                cutoff <- min(head(abs_logr, n=n))
            }
            logr2 <- logr[abs(logr) < cutoff]
            v <- unname(as.numeric(var(logr2)))
            list(
            mean     = unname(as.numeric(mean(logr2))),
            stdev    = sqrt(v),
            var      = v,
            skewness = unname(as.numeric(moments::skewness(logr2))),
            kurtosis = unname(as.numeric(moments::kurtosis(logr2)))
            )
        }
        log2n_wide <- log(n)/log(2)
        asymp_stats_wide <- lapply(n-1,stats_n)
        asymp_stdev_wide <- sapply(n, function(i) asymp_stats_wide[[i]]$stdev)
        asymp_skewness_wide <- sapply(n, function(i) asymp_stats_wide[[i]]$skewness)
        asymp_kurt_wide <- sapply(n, function(i) asymp_stats_wide[[i]]$kurtosis)

        i <- 2:(2^tail.N2)
        log2n <- log2n_wide[i]
        asymp_stdev <- asymp_stdev_wide[i]
        fit_asymp_stdev <- lm(asymp_stdev ~ poly(log2n,2,raw=TRUE))
        asymp_skewness <- asymp_skewness_wide[i]
        fit_asymp_skewness <- lm(asymp_skewness ~ log2n)
        asymp_kurt <- asymp_kurt_wide[i]
        fit_asymp_kurt <- lm(asymp_kurt ~ log2n)
        
        log2n_line <- seq(0, tail.N2, by=0.1)
        asymp_stdev_line <- predict(fit_asymp_stdev, data.frame(log2n=c(log2n_line)))
        asymp_skewness_line <- predict(fit_asymp_skewness, data.frame(log2n=c(log2n_line)))
        asymp_kurt_line <- predict(fit_asymp_kurt, data.frame(log2n=c(log2n_line)))
        
        asymp_stdev_0 <- predict(fit_asymp_stdev, data.frame(log2n=c(0)))
        asymp_stdev_1 <- predict(fit_asymp_stdev, data.frame(log2n=c(1)))
        asymp_skewness_0 <- predict(fit_asymp_skewness, data.frame(log2n=c(0)))
        asymp_skewness_1 <- predict(fit_asymp_skewness, data.frame(log2n=c(1)))
        asymp_kurt_0 <- predict(fit_asymp_kurt, data.frame(log2n=c(0)))
        asymp_kurt_1 <- predict(fit_asymp_kurt, data.frame(log2n=c(1)))

        tail <- list(
            N1 = tail.N1,
            N2 = tail.N2,
            
            log2n_wide = log2n_wide,
            asymp_stdev_wide = asymp_stdev_wide,
            asymp_skewness_wide = asymp_skewness_wide,
            asymp_kurt_wide = asymp_kurt_wide,
            
            log2n = log2n,
            asymp_stdev = asymp_stdev,
            fit_asymp_stdev = fit_asymp_stdev,
            asymp_skewness = asymp_skewness,
            fit_asymp_skewness = fit_asymp_skewness,
            asymp_kurt = asymp_kurt,
            fit_asymp_kurt = fit_asymp_kurt,
            
            log2n_line = log2n_line,
            asymp_stdev_line = asymp_stdev_line,
            asymp_skewness_line = asymp_skewness_line,
            asymp_kurt_line = asymp_kurt_line,
            asymp_stdev_0 = asymp_stdev_0,
            asymp_stdev_1 = asymp_stdev_1,
            asymp_skewness_0 = asymp_skewness_0,
            asymp_skewness_1 = asymp_skewness_1,
            asymp_kurt_0 = asymp_kurt_0,
            asymp_kurt_1 = asymp_kurt_1
        )
        xtsAttributes(ts) <- list(tail = tail)
    }
    ts
}

