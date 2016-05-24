#' Conditionally detrend a time series
#' 
#' This function detrends a time series when its linear trend is more significant than a threshold.
#' 
#' @param x a vector or time series.
#' @param thres a scalar specifying the threshold. When the adjusted R square coefficient of the linear fitting 
#' is larger than this threshold, the linear trend is substracted from the original time series. Default is 0.85.
#' @return detrended \code{x}.
#' @export
#' @examples
#' t=seq(0.001,1,0.001)
#' set.seed(123)
#' x=10*t+rnorm(1000)
#' dtrx=detrendc(x)
#' # plot the simulated x
#' plot(t,x,ty='l',xlab='t',ylab='x')
#' # plot the detrended x
#' lines(t,dtrx,col=2)
#' legend(0,12,legend=c('x','detrended x'),col=c(1,2),lty=1)


detrendc <- function(x, thres = 0.85) {
    n = length(x)
    t = seq(1, n)
    # fit linear regression
    model = summary(fastLm(cbind(1, t), x))
    coef = model$coefficients[2, 1]
    # if the linear trend is more significant than the threshold, it is substracted from x
    if (model$adj.r.squared > thres) {
        detrended_x = x - coef * t
    } else {
        detrended_x = x
    }
    return(detrended_x)
} 
