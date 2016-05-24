#' Description: Sensitivity Early Warning Signals
#'
#' \code{sensitivity_ews} is used to estimate trends in statistical moments for different sizes of rolling windows along a timeseries. The trends are estimated by the nonparametric Kendall tau correlation coefficient.
#'
# Details:
#' see ref below
#'
#' Arguments:
#'    @param timeseries a numeric vector of the observed univariate timeseries values or a numeric matrix where the first column represents the time index and the second the observed timeseries values. Use vectors/matrices with headings.
#'    @param indicator is the statistic (leading indicator) selected for which the sensitivity analysis is perfomed. Currently, the indicators supported are: \code{ar1} autoregressive coefficient of a first order AR model, \code{sd} standard deviation, \code{acf1} autocorrelation at first lag, \code{sk} skewness, \code{kurt} kurtosis, \code{cv} coeffcient of variation, \code{returnrate}, and \code{densratio} density ratio of the power spectrum at low frequencies over high frequencies.
#'    @param winsizerange is the range of the rolling window sizes expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default is 25\% - 75\%.
#'    @param incrwinsize increments the rolling window size (must be numeric between 0 and 100). Default is 25.
#'    @param detrending the timeseries can be detrended/filtered. There are three options: \code{gaussian} filtering, \code{loess} fitting, \code{linear} detrending and \code{first-differencing}. Default is \code{no} detrending.
#'    @param bandwidthrange is the range of the bandwidth used for the Gaussian kernel when gaussian filtering is selected. It is expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default is 5\% - 100\%.
#'    @param spanrange parameter that controls the degree of smoothing (numeric between 0 and 100). Default is 5\% - 100\%. see more on loess{stats}
#'    @param degree the degree of polynomial to be used for when loess fitting is applied, normally 1 or 2 (Default). see more on loess{stats}
#'    @param incrbandwidth is the size to increment the bandwidth used for the Gaussian kernel when gaussian filtering is applied. It is expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default is 20.
#'    @param incrspanrange Span range
#'    @param logtransform logical. If TRUE data are logtransformed prior to analysis as log(X+1). Default is FALSE. 
#'    @param interpolate logical. If TRUE linear interpolation is applied to produce a timeseries of equal length as the original. Default is FALSE (assumes there are no gaps in the timeseries).
#' 
# Returns:
#'  @return \code{sensitivity_ews} returns a matrix that contains the Kendall tau rank correlation estimates for the rolling window sizes (rows) and bandwidths (columns), if \code{gaussian filtering} is selected. 
#'  
#'  In addition, \code{sensitivity_ews} returns a plot with the Kendall tau estimates and their p-values for the range of rolling window sizes used, together with a histogram of the distributions of the statistic and its significance. When \code{gaussian filtering} is chosen, a contour plot is produced for the Kendall tau estimates and their p-values for the range of both rolling window sizes and bandwidth used. A reverse triangle indicates the combination of the two parameters for which the Kendall tau was the highest
#'  
#' @export
#' 
#' @importFrom moments skewness
#' @importFrom moments kurtosis
#' @importFrom fields image.plot
#'
#' @author Vasilis Dakos \email{vasilis.dakos@@gmail.com}
#' @references Dakos, V., et al (2008). 'Slowing down as an early warning signal for abrupt climate change.' \emph{Proceedings of the National Academy of Sciences} 105(38): 14308-14312 
#' 
#' Dakos, V., et al (2012).'Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data.' \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010 
#' @seealso 
#' \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}}; \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}}; \code{\link{surrogates_ews}}; \code{\link{ch_ews}}; \code{\link{movpotential_ews}}; \code{\link{livpotential_ews}}
# ; \code{\link{timeVAR_ews}}; \code{\link{thresholdAR_ews}}
#' @examples 
#' data(foldbif)
#' output=sensitivity_ews(foldbif,indicator='sd',detrending='gaussian',
#' incrwinsize=25,incrbandwidth=20)
#' @keywords early-warning

# Author: Vasilis Dakos, January 4, 2012

sensitivity_ews <- function(timeseries, indicator = c("ar1", "sd", "acf1", "sk", 
    "kurt", "cv", "returnrate", "densratio"), winsizerange = c(25, 75), incrwinsize = 25, 
    detrending = c("no", "gaussian", "loess", "linear", "first-diff"), bandwidthrange = c(5, 
        100), spanrange = c(5, 100), degree = NULL, incrbandwidth = 20, incrspanrange = 10, 
    logtransform = FALSE, interpolate = FALSE) {
    
    # timeseries<-ts(timeseries) #strict data-types the input data as tseries object
    # for use in later steps
    
    timeseries <- data.matrix(timeseries)
    if (dim(timeseries)[2] == 1) {
        Y = timeseries
        timeindex = 1:dim(timeseries)[1]
    } else if (dim(timeseries)[2] == 2) {
        Y <- timeseries[, 2]
        timeindex <- timeseries[, 1]
    } else {
        warning("not right format of timeseries input")
    }
    
    # Interpolation
    if (interpolate) {
        YY <- approx(timeindex, Y, n = length(Y), method = "linear")
        Y <- YY$y
    } else {
        Y <- Y
    }
    
    # Log-transformation
    if (logtransform) {
        Y <- log(Y + 1)
    }
    
    # Determine the step increases in rolling windowsize
    incrtw <- incrwinsize
    tw <- seq(floor(winsizerange[1] * length(Y)/100), floor(winsizerange[2] * length(Y)/100), 
        by = incrtw)
    twcol <- length(tw)
    low <- 6
    
    # Detrending
    detrending <- match.arg(detrending)
    if (detrending == "gaussian") {
        incrbw <- incrbandwidth
        width <- seq(floor(bandwidthrange[1] * length(Y)/100), floor(bandwidthrange[2] * 
            length(Y)/100), by = incrbw)
        bwrow <- length(width)
        # Create matrix to store Kendall trend statistics
        Ktauestind <- matrix(, bwrow, twcol)
        Ktaupind <- matrix(, bwrow, twcol)
        # Estimation
        for (wi in 1:(length(width))) {
            width1 <- width[wi]
            smYY <- ksmooth(timeindex, Y, kernel = c("normal"), bandwidth = width1, 
                range.x = range(timeindex), n.points = length(timeindex))
            nsmY <- Y - smYY$y
            for (ti in 1:length(tw)) {
                tw1 <- tw[ti]
                # Rearrange data for indicator calculation
                omw1 <- length(nsmY) - tw1 + 1  ##number of overlapping moving windows
                high <- omw1
                nMR1 <- matrix(data = NA, nrow = tw1, ncol = omw1)
                for (i in 1:omw1) {
                  Ytw <- nsmY[i:(i + tw1 - 1)]
                  nMR1[, i] <- Ytw
                }
                # Estimate indicator
                indicator = match.arg(indicator)
                if (indicator == "ar1") {
                  indic <- apply(nMR1, 2, function(x) {
                    nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, dmean = FALSE, 
                      intercept = FALSE)
                    nAR1$ar
                  })
                } else if (indicator == "sd") {
                  indic <- apply(nMR1, 2, sd)
                } else if (indicator == "sk") {
                  indic <- apply(nMR1, 2, skewness)
                } else if (indicator == "kurt") {
                  indic <- apply(nMR1, 2, kurtosis)
                } else if (indicator == "acf1") {
                  indic <- apply(nMR1, 2, function(x) {
                    nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
                    nACF$acf[2]
                  })
                } else if (indicator == "returnrate") {
                  indic <- apply(nMR1, 2, function(x) {
                    nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
                    1 - nACF$acf[2]
                  })
                } else if (indicator == "cv") {
                  indic <- apply(nMR1, 2, function(x) {
                    sd(x)/mean(x)
                  })
                } else if (indicator == "densratio") {
                  indic <- apply(nMR1, 2, function(x) {
                    spectfft <- spec.ar(x, n.freq = omw1, plot = FALSE, order = 1)
                    spectfft$spec
                    spectfft$spec[low]/spectfft$spec[high]
                  })
                }
                # Calculate trend statistics
                timevec <- seq(1, length(indic))
                Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
                  conf.level = 0.95)
                Ktauestind[wi, ti] <- Kt$estimate
                Ktaupind[wi, ti] <- Kt$p.value
            }
        }
        # Plot
        layout(matrix(1:4, 2, 2))
        par(font.main = 10, mar = (c(4.6, 4, 0.5, 2) + 0.2), mgp = c(2, 1, 0), oma = c(0.5, 
            1, 2, 1))
        image.plot(width, tw, Ktauestind, zlim = c(-1, 1), xlab = "filtering bandwidth", 
            ylab = "rolling window size", main = "Kendall tau estimate", cex.main = 0.8, 
            log = "y", nlevel = 20, col = rainbow(20))
        ind = which(Ktauestind == max(Ktauestind), arr.ind = TRUE)
        lines(width[ind[1]], tw[ind[2]], type = "p", cex = 1.2, pch = 17, col = 1)
        hist(Ktauestind, breaks = 12, col = "lightblue", main = NULL, xlab = "Kendall tau estimate", 
            ylab = "occurence", border = "black", xlim = c(-1, 1))
        fields::image.plot(width, tw, Ktaupind, zlim = c(0, max(Ktaupind)), xlab = "filtering bandwidth", 
            ylab = "rolling window size", main = "Kendall tau p-value", log = "y", 
            cex.main = 0.8, nlevel = 20, col = rainbow(20))
        lines(width[ind[1]], tw[ind[2]], type = "p", cex = 1.2, pch = 17, col = 1)
        hist(Ktaupind, breaks = 12, col = "yellow", main = NULL, xlab = "Kendall tau p-value", 
            ylab = "occurence", border = "black", xlim = c(0, max(Ktaupind)))
        mtext(paste("Indicator ", toupper(indicator)), side = 3, line = 0.2, outer = TRUE)
        
        # Output
        out <- data.frame(Ktauestind)
        colnames(out) <- tw
        rownames(out) <- width
        return(out)
        
    } else if (detrending == "loess") {
        incrbw <- incrspanrange
        width <- seq(floor(spanrange[1] * length(Y)/100), floor(spanrange[2] * length(Y)/100), 
            by = incrbw)
        bwrow <- length(width)
        # Create matrix to store Kendall trend statistics
        Ktauestind <- matrix(, bwrow, twcol)
        Ktaupind <- matrix(, bwrow, twcol)
        # Estimation
        if (is.null(degree)) {
            degree <- 2
        } else {
            degree <- degree
        }
        for (wi in 1:(length(width))) {
            width1 <- width[wi]
            smYY <- loess(Y ~ timeindex, span = width1, degree = degree, normalize = FALSE, 
                family = "gaussian")
            smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
            nsmY <- Y - smY
            for (ti in 1:length(tw)) {
                tw1 <- tw[ti]
                # Rearrange data for indicator calculation
                omw1 <- length(nsmY) - tw1 + 1  ##number of overlapping moving windows
                high <- omw1
                nMR1 <- matrix(data = NA, nrow = tw1, ncol = omw1)
                for (i in 1:omw1) {
                  Ytw <- nsmY[i:(i + tw1 - 1)]
                  nMR1[, i] <- Ytw
                }
                # Estimate indicator
                indicator = match.arg(indicator)
                if (indicator == "ar1") {
                  indic <- apply(nMR1, 2, function(x) {
                    nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, dmean = FALSE, 
                      intercept = FALSE)
                    nAR1$ar
                  })
                } else if (indicator == "sd") {
                  indic <- apply(nMR1, 2, sd)
                } else if (indicator == "sk") {
                  indic <- apply(nMR1, 2, skewness)
                } else if (indicator == "kurt") {
                  indic <- apply(nMR1, 2, kurtosis)
                } else if (indicator == "acf1") {
                  indic <- apply(nMR1, 2, function(x) {
                    nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
                    nACF$acf[2]
                  })
                } else if (indicator == "returnrate") {
                  indic <- apply(nMR1, 2, function(x) {
                    nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
                    1 - nACF$acf[2]
                  })
                } else if (indicator == "cv") {
                  indic <- apply(nMR1, 2, function(x) {
                    sd(x)/mean(x)
                  })
                } else if (indicator == "densratio") {
                  indic <- apply(nMR1, 2, function(x) {
                    spectfft <- spec.ar(x, n.freq = omw1, plot = FALSE, order = 1)
                    spectfft$spec
                    spectfft$spec[low]/spectfft$spec[high]
                  })
                }
                # Calculate trend statistics
                timevec <- seq(1, length(indic))
                Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
                  conf.level = 0.95)
                Ktauestind[wi, ti] <- Kt$estimate
                Ktaupind[wi, ti] <- Kt$p.value
            }
        }
        # Plot
        layout(matrix(1:4, 2, 2))
        par(font.main = 10, mar = (c(4.6, 4, 0.5, 2) + 0.2), mgp = c(2, 1, 0), oma = c(0.5, 
            1, 2, 1))
        fields::image.plot(width, tw, Ktauestind, zlim = c(-1, 1), xlab = "filtering bandwidth", 
            ylab = "rolling window size", main = "Kendall tau estimate", cex.main = 0.8, 
            log = "y", nlevel = 20, col = rainbow(20))
        ind = which(Ktauestind == max(Ktauestind), arr.ind = TRUE)
        lines(width[ind[1]], tw[ind[2]], type = "p", cex = 1.2, pch = 17, col = 1)
        hist(Ktauestind, breaks = 12, col = "lightblue", main = NULL, xlab = "Kendall tau estimate", 
            ylab = "occurence", border = "black", xlim = c(-1, 1))
        fields::image.plot(width, tw, Ktaupind, zlim = c(0, max(Ktaupind)), xlab = "filtering bandwidth", 
            ylab = "rolling window size", main = "Kendall tau p-value", log = "y", 
            cex.main = 0.8, nlevel = 20, col = rainbow(20))
        lines(width[ind[1]], tw[ind[2]], type = "p", cex = 1.2, pch = 17, col = 1)
        hist(Ktaupind, breaks = 12, col = "yellow", main = NULL, xlab = "Kendall tau p-value", 
            ylab = "occurence", border = "black", xlim = c(0, max(Ktaupind)))
        mtext(paste("Indicator ", toupper(indicator)), side = 3, line = 0.2, outer = TRUE)
        
        # Output
        out <- data.frame(Ktauestind)
        colnames(out) <- tw
        rownames(out) <- width
        return(out)
        
    } else if (detrending == "linear") {
        nsmY <- resid(lm(Y ~ timeindex))
    } else if (detrending == "first-diff") {
        nsmY <- diff(Y)
    } else if (detrending == "no") {
        nsmY <- Y
    }
    
    # Create matrix to store Kendall trend statistics
    Ktauestind <- matrix(, twcol, 1)
    Ktaupind <- matrix(, twcol, 1)
    
    for (ti in 1:length(tw)) {
        tw1 <- tw[ti]
        # Rearrange data for indicator calculation
        omw1 <- length(nsmY) - tw1 + 1  ##number of overlapping moving windows
        high = omw1
        nMR1 <- matrix(data = NA, nrow = tw1, ncol = omw1)
        for (i in 1:omw1) {
            Ytw <- nsmY[i:(i + tw1 - 1)]
            nMR1[, i] <- Ytw
        }
        # Estimate indicator
        indicator = match.arg(indicator)
        if (indicator == "ar1") {
            indic <- apply(nMR1, 2, function(x) {
                nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)
                nAR1$ar
            })
        } else if (indicator == "sd") {
            indic <- apply(nMR1, 2, sd)
        } else if (indicator == "sk") {
            indic <- apply(nMR1, 2, skewness)
        } else if (indicator == "kurt") {
            indic <- apply(nMR1, 2, kurtosis)
        } else if (indicator == "acf1") {
            indic <- apply(nMR1, 2, function(x) {
                nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
                nACF$acf[2]
            })
        } else if (indicator == "returnrate") {
            indic <- apply(nMR1, 2, function(x) {
                nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
                1 - nACF$acf[2]
            })
        } else if (indicator == "cv") {
            indic <- apply(nMR1, 2, function(x) {
                sd(x)/mean(x)
            })
        } else if (indicator == "densratio") {
            indic <- apply(nMR1, 2, function(x) {
                spectfft <- spec.ar(x, n.freq = omw1, plot = FALSE, order = 1)
                spectfft$spec
                spectfft$spec[low]/spectfft$spec[high]
            })
        }
        # Calculate trend statistics
        timevec <- seq(1, length(indic))
        Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
            conf.level = 0.95)
        Ktauestind[ti] <- Kt$estimate
        Ktaupind[ti] <- Kt$p.value
    }
    # Plotting
    dev.new()
    layout(matrix(1:4, 2, 2))
    par(font.main = 10, mar = (c(4.6, 4, 0.5, 2) + 0.2), mgp = c(2, 1, 0), oma = c(0.5, 
        1, 2, 1))
    plot(tw, Ktauestind, type = "l", ylim = c(-1, 1), log = "x", xlab = "rolling window size", 
        ylab = "Kendall tau estimate")
    hist(Ktauestind, breaks = 12, col = "lightblue", xlab = "Kendall tau estimate", 
        ylab = "occurence", border = "black", xlim = c(-1, 1), main = NULL)
    plot(tw, Ktaupind, type = "l", xlab = "rolling window size", log = "x", ylab = "Kendall tau p-value", 
        ylim = c(0, max(Ktaupind)))
    hist(Ktaupind, breaks = 12, col = "yellow", xlab = "Kendall tau p-value", main = NULL, 
        ylab = "occurence", border = "black", xlim = c(0, max(Ktaupind)))
    mtext(paste("Indicator ", toupper(indicator)), side = 3, line = 0.2, outer = TRUE)
    
    # Output out<-data.frame(tw,Ktauestind,Ktaupind) colnames(out)<-c('rolling
    # window','Kendall tau estimate','Kendall tau p-value')
    out <- data.frame(Ktauestind)
    rownames(out) <- tw
    return(out)
} 
