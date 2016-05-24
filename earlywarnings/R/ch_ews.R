#' Description: Conditional Heteroskedasticity
#'
#' \code{ch_ews} is used to estimate changes in conditional heteroskedasticity within rolling windows along a timeseries
#'
# Details:
#' see ref below
#'
#' Arguments:
#'    @param timeseries a numeric vector of the observed timeseries values or a numeric matrix where the first column represents the time index and the second the observed timeseries values. Use vectors/matrices with headings.
#'    @param winsize is length of the rolling window expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default is 10\%.
#'    @param alpha is the significance threshold (must be numeric). Default is 0.1.
#'    @param optim logical. If TRUE an autoregressive model is fit to the data within the rolling window using AIC optimization. Otherwise an autoregressive model of specific order \code{lags} is selected.
#'    @param lags is a parameter that determines the specific order of an autoregressive model to fit the data. Default is 4.
#'    @param logtransform logical. If TRUE data are logtransformed prior to analysis as log(X+1). Default is FALSE.
#'    @param interpolate logical. If TRUE linear interpolation is applied to produce a timeseries of equal length as the original. Default is FALSE (assumes there are no gaps in the timeseries). 
#' 
# Returns:
#'   @return \code{ch_ews} returns a matrix that contains:
#'   @return \item{time}{the time index.}
#'   @return \item{r.squared}{the R2 values of the regressed residuals.}
#'   @return \item{critical.value}{the chi-square critical value based on the desired \code{alpha} level for 1 degree of freedom divided by the number of residuals used in the regression.}
#'   @return \item{test.result}{logical. It indicates whether conditional heteroskedasticity was significant.}
#'   @return \item{ar.fit.order}{the order of the specified autoregressive model- only informative if \code{optim} FALSE was selected.}
#'
#' In addition, \code{ch_ews} plots the original timeseries and the R2 where the level of significance is also indicated.
#'  
#' @export
#' 
#' @author T. Cline, modified by V. Dakos
#' @references Seekell, D. A., et al (2011). 'Conditional heteroscedasticity as a leading indicator of ecological regime shifts.' \emph{American Naturalist} 178(4): 442-451
#' 
#' Dakos, V., et al (2012).'Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data.' \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010 
#' @seealso 
#' \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}}; \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}}; \code{\link{surrogates_ews}}; \code{\link{ch_ews}}; \code{movpotential_ews}; \code{livpotential_ews} 
# ; \code{\link{timeVAR_ews}}; \code{\link{thresholdAR_ews}}
#' @examples 
#' data(foldbif)
#' out=ch_ews(foldbif, winsize=50, alpha=0.05, optim=TRUE, lags)
#' @keywords early-warning
# Author: Timothy Cline, October 25, 2011.  Modified by: Vasilis Dakos, January
# 3, 2012.

ch_ews <- function(timeseries, winsize = 10, alpha = 0.1, optim = TRUE, lags = 4, 
    logtransform = FALSE, interpolate = FALSE) {
    
    timeseries <- ts(timeseries)  #strict data-types the input data as tseries object for use in later steps
    if (dim(timeseries)[2] == 1) {
        ts.in = timeseries
        timeindex = 1:dim(timeseries)[1]
    } else if (dim(timeseries)[2] == 2) {
        ts.in = timeseries[, 2]
        timeindex = timeseries[, 1]
    } else {
        warning("not right format of timeseries input")
    }
    
    # Interpolation
    if (interpolate) {
        YY <- approx(timeindex, ts.in, n = length(ts.in), method = "linear")
        ts.in <- YY$y
    }
    
    # Log-transformation
    if (logtransform) {
        ts.in <- log(ts.in + 1)
    }
    
    winSize = round(winsize * length(ts.in)/100)
    sto <- matrix(nrow = (length(ts.in) - (winSize - 1)), ncol = 5)  # creates a matrix to store output
    
    count <- 1  #place holder for writing to the matrix
    for (i2 in winSize:length(ts.in)) {
        # loop to iterate through the model values by window lengths of the input value
        
        # the next line applys the autoregressive model optimized using AIC then we omit
        # the first data point(s) which come back as NA and square the residuals
        if (optim == TRUE) {
            arm <- ar(ts.in[(i2 - (winSize - 1)):i2], method = "ols")
        } else {
            arm <- ar(ts.in[(i2 - (winSize - 1)):i2], aic = FALSE, order.max = lags, 
                method = "ols")
        }
        resid1 <- na.omit(arm$resid)^2
        
        l1 <- length(resid1)  # stores the number of residuals for many uses later
        lm1 <- lm(resid1[2:l1] ~ resid1[1:(l1 - 1)])  #calculates simple OLS model of describing t+1 by t
        
        # calculates the critical value: Chi-squared critical value using desired alpha
        # level and 1 degree of freedom / number of residuals used in regression
        critical <- qchisq((1 - alpha), df = 1)/(length(resid1) - 1)
        
        sto[count, 1] <- timeindex[i2]  # stores a time component
        sto[count, 2] <- summary(lm1)$r.squared  # stores the r.squared for this window
        sto[count, 3] <- critical  # stores the critical value
        
        # the next flow control group stores a simple 1 for significant test or 0 for
        # non-significant test
        if (summary(lm1)$r.squared > critical) {
            sto[count, 4] <- 1
        } else {
            sto[count, 4] <- 0
        }
        sto[count, 5] <- arm$order
        count <- count + 1  # increment the place holder
    }
    
    sto <- data.frame(sto)  # data types the matrix as a data frame
    colnames(sto) <- c("time", "r.squared", "critical.value", "test.result", "ar.fit.order")  # applies column names to the data frame
    
    # This next series of flow control statements will adjust the max and minimum
    # values to yield prettier plots In some circumstances it is possible for all
    # values to be far greater or far less than the critical value; in all cases we
    # want the critical line ploted on the figure
    if (max(sto$r.squared) < critical) {
        maxY <- critical + 0.02
        minY <- critical - 0.02
    } else if (min(sto$r.squared >= critical)) {
        minY <- critical - 0.02
        maxY <- critical + 0.02
    } else {
        maxY <- max(sto$r.squared)
        minY <- min(sto$r.squared)
    }
    
    # this creates a very simple plot that is well fitted to the data.  it also plots
    # the critical value line
    
    par(mar = (c(0, 4, 0, 1) + 0), oma = c(5, 1, 2, 1), mfrow = c(2, 1))
    plot(timeindex, ts.in, type = "l", ylab = "data", xlab = "", cex.axis = 0.8, 
        cex.lab = 0.8, xaxt = "n", las = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))
    plot(timeindex[winSize:length(timeindex)], sto$r.squared, ylab = expression(paste("R^2")), 
        xlab = "time", type = "b", cex = 0.5, cex.lab = 0.8, cex.axis = 0.8, las = 1, 
        ylim = c(min(sto$r.squared), max(sto$r.squared)), xlim = c(timeindex[1], 
            timeindex[length(timeindex)]))
    legend("topleft", "conditional heteroskedasticity", bty = "n")
    abline(h = sto$critical, lwd = 0.5, lty = 2, col = 2)
    mtext("time", side = 1, line = 2, cex = 0.8)  #outer=TRUE print on the outer margin
    
    return(sto)
} 
