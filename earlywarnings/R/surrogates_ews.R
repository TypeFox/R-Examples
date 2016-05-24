#' Description: Surrogates Early Warning Signals
#'
#' \code{surrogates_ews} is used to estimate distributions of trends in statistical moments from different surrogate timeseries generated after fitting an ARMA(p,q) model on the data. The trends are estimated by the nonparametric Kendall tau correlation coefficient and can be compared to the trends estimated in the original timeseries to produce probabilities of false positives.
#'
# Details:
#' see ref below
#'
#' Arguments:
#'    @param timeseries a numeric vector of the observed univariate timeseries values or a numeric matrix where the first column represents the time index and the second the observed timeseries values. Use vectors/matrices with headings.
#'    @param indicator is the statistic (leading indicator) selected for which the surrogate timeseries are produced. Currently, the indicators supported are: \code{ar1} autoregressive coefficient of a first order AR model, \code{sd} standard deviation, \code{acf1} autocorrelation at first lag, \code{sk} skewness, \code{kurt} kurtosis, \code{cv} coeffcient of variation, \code{returnrate}, and \code{densratio} density ratio of the power spectrum at low frequencies over high frequencies.
#'    @param winsize is the size of the rolling window expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default valuise 50\%.
#'    @param detrending the timeseries can be detrended/filtered prior to analysis. There are three options: \code{gaussian} filtering, \code{loess} fitting, \code{linear} detrending and \code{first-diff}erencing. Default is \code{no} detrending.
#'    @param bandwidth is the bandwidth used for the Gaussian kernel when gaussian filtering is selected. It is expressed as percentage of the timeseries length (must be numeric between 0 and 100). Alternatively it can be given by the bandwidth selector \code{\link{bw.nrd0}} (Default).
#'    @param span parameter that controls the degree of smoothing (numeric between 0 and 100, Default 25). see more on loess{stats}
#'    @param degree the degree of polynomial to be used for when loess fitting is applied, normally 1 or 2 (Default). see more on loess{stats}
#'    @param boots the number of surrogate data. Default is 100.
#'    @param logtransform logical. If TRUE data are logtransformed prior to analysis as log(X+1). Default is FALSE.
#'    @param interpolate logical. If TRUE linear interpolation is applied to produce a timeseries of equal length as the original. Default is FALSE (assumes there are no gaps in the timeseries). 
#' 
# Returns:
#'   @return \code{surrogates_ews} returns a matrix that contains:
#'   @return \item{Kendall tau estimate original}{the trends of the original timeseries.}
#'   @return \item{Kendall tau p-value original}{the p-values of the trends of the original timeseries.}
#'   @return \item{Kendall tau estimate surrogates}{the trends of the surrogate timeseries.}
#'   @return \item{Kendall tau p-value surrogates}{the associated p-values of the trends of the surrogate timeseries.}
#'   @return \item{significance p}{the p-value for the original Kendall tau rank correlation estimate compared to the surrogates.}
#'
#' In addition, \code{surrogates_ews} returns a plot with the distribution of the surrogate Kendall tau estimates and the Kendall tau estimate of the original series. Vertical lines indicate the 5\% and 95\% significance levels.
#'  
#' @export
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
#' output=surrogates_ews(foldbif,indicator='sd',winsize=50,detrending='gaussian',
#' bandwidth=10,boots=200,logtransform=FALSE,interpolate=FALSE)
#' @keywords early-warning

# Author: Vasilis Dakos, January 4, 2012

surrogates_ews <- function(timeseries, indicator = c("ar1", "sd", "acf1", "sk", "kurt", 
    "cv", "returnrate", "densratio"), winsize = 50, detrending = c("no", "gaussian", 
    "loess", "linear", "first-diff"), bandwidth = NULL, span = NULL, degree = NULL, 
    boots = 100, logtransform = FALSE, interpolate = FALSE) {
    
    # timeseries<-ts(timeseries) #strict data-types the input data as tseries object
    # for use in later steps
    
    skewness <- moments::skewness
    kurtosis <- moments::kurtosis
    
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
    
    # Detrending
    detrending <- match.arg(detrending)
    if (detrending == "gaussian") {
        if (is.null(bandwidth)) {
            bw <- round(bw.nrd0(timeindex))
        } else {
            bw <- round(length(Y) * bandwidth)/100
        }
        smYY <- ksmooth(timeindex, Y, kernel = c("normal"), bandwidth = bw, range.x = range(timeindex), 
            n.points = length(timeindex))
        nsmY <- Y - smYY$y
        smY <- smYY$y
    } else if (detrending == "loess") {
        if (is.null(span)) {
            span <- 25/100
        } else {
            span <- span/100
        }
        if (is.null(degree)) {
            degree <- 2
        } else {
            degree <- degree
        }
        smYY <- loess(Y ~ timeindex, span = span, degree = degree, normalize = FALSE, 
            family = "gaussian")
        smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
        nsmY <- Y - smY
    } else if (detrending == "linear") {
        nsmY <- resid(lm(Y ~ timeindex))
        smY <- fitted(lm(Y ~ timeindex))
    } else if (detrending == "first-diff") {
        nsmY <- diff(Y)
        timeindexdiff <- timeindex[1:(length(timeindex) - 1)]
    } else if (detrending == "no") {
        smY <- Y
        nsmY <- Y
    }
    
    
    # Rearrange data for indicator calculation
    mw <- round(length(Y) * winsize)/100
    omw <- length(nsmY) - mw + 1
    low <- 6
    high <- omw
    nMR <- matrix(data = NA, nrow = mw, ncol = omw)
    for (i in 1:omw) {
        Ytw <- nsmY[i:(i + mw - 1)]
        nMR[, i] <- Ytw
    }
    # Estimate indicator
    indicator = match.arg(indicator)
    if (indicator == "ar1") {
        indic <- apply(nMR, 2, function(x) {
            nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)
            nAR1$ar
        })
    } else if (indicator == "sd") {
        indic <- apply(nMR, 2, sd)
    } else if (indicator == "sk") {
        indic <- apply(nMR, 2, skewness)
    } else if (indicator == "kurt") {
        indic <- apply(nMR, 2, kurtosis)
    } else if (indicator == "acf1") {
        indic <- apply(nMR, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            nACF$acf[2]
        })
    } else if (indicator == "returnrate") {
        indic <- apply(nMR, 2, function(x) {
            nACF <- acf(x, lag.max = 1, type = c("correlation"), plot = FALSE)
            1 - nACF$acf[2]
        })
    } else if (indicator == "cv") {
        indic <- apply(nMR, 2, function(x) {
            sd(x)/mean(x)
        })
    } else if (indicator == "densratio") {
        indic <- apply(nMR, 2, function(x) {
            spectfft <- spec.ar(x, n.freq = omw, plot = FALSE, order = 1)
            spectfft$spec
            spectfft$spec[low]/spectfft$spec[high]
        })
    }
    # Calculate trend statistics
    timevec <- seq(1, length(indic))
    Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
        conf.level = 0.95)
    Ktauestindorig <- Kt$estimate
    Ktaupindorig <- Kt$p.value
    
    # Fit ARMA model based on AIC
    arma = matrix(, 4, 5)
    for (ij in 1:4) {
        for (jj in 0:4) {
            ARMA <- arima(nsmY, order = c(ij, 0, jj), include.mean = FALSE)
            arma[ij, jj + 1] = ARMA$aic
            print(paste("AR", "MA", "AIC"), quote = FALSE)
            print(paste(ij, jj, ARMA$aic), zero.print = ".", quote = FALSE)
        }
    }
    
    # Simulate ARMA(p,q) model fitted on residuals
    ind = which(arma == min(arma), arr.ind = TRUE)
    ARMA <- arima(nsmY, order = c(ind[1], 0, ind[2] - 1), include.mean = FALSE)
    
    Ktauestind <- numeric()
    Ktaupind <- numeric()
    
    for (jjj in 1:boots) {
        x = arima.sim(n = length(nsmY), list(ar = c(ARMA$coef[1:ind[1]]), ma = c(ARMA$coef[(1 + 
            ind[1]):(ind[1] + ind[2] - 1)])), sd = sqrt(ARMA$sigma2))
        
        ## Rearrange data for indicator calculation
        nMR1 <- matrix(data = NA, nrow = mw, ncol = omw)
        for (i in 1:omw) {
            Ytw <- x[i:(i + mw - 1)]
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
                spectfft <- spec.ar(x, n.freq = omw, plot = FALSE, order = 1)
                spectfft$spec
                spectfft$spec[low]/spectfft$spec[high]
            })
        }
        
        # Calculate trend statistics
        timevec <- seq(1, length(indic))
        Kt <- cor.test(timevec, indic, alternative = c("two.sided"), method = c("kendall"), 
            conf.level = 0.95)
        Ktauestind[jjj] <- Kt$estimate
        Ktaupind[jjj] <- Kt$p.value
    }
    
    # Estimate probability of false positive
    q <- sort(Ktauestind, na.last = NA)
    Kpos <- max(which(Ktauestindorig > q), na.rm = TRUE)
    p <- (boots + 1 - Kpos)/boots
    print(paste("significance p = ", p, " estimated from ", boots, " surrogate ARMA timeseries"))
    
    # Plotting
    dev.new()
    par(font.main = 10, mar = (c(4.6, 4, 0.5, 2) + 0.2), oma = c(0.5, 1, 2, 1))
    hist(Ktauestind, freq = TRUE, nclass = 20, xlim = c(-1, 1), col = "blue", main = NULL, 
        xlab = "Surrogate Kendall tau estimates", ylab = "occurrence", ylim = c(0, 
            boots))
    abline(v = q[0.05 * boots], col = "black", lwd = 1)
    abline(v = q[0.95 * boots], col = "black", lwd = 1)
    points(Ktauestindorig, 0, pch = 21, bg = "black", col = "black", cex = 1)
    mtext(paste("Indicator ", toupper(indicator)), side = 3, line = 0.2, outer = TRUE)
    
    # Output
    out <- data.frame(Ktauestindorig, Ktaupindorig, Ktauestind, Ktaupind, p)
    colnames(out) <- c("Kendall tau estimate original", "Kendall tau p-value original", 
        "Kendall tau estimate surrogates", "Kendall tau p-value surrogates", "significance p")
    return(out)
} 
