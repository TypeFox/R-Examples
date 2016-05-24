#' Description: Generic Early Warning Signals
#'
#' \code{generic_ews} is used to estimate statistical moments within rolling windows along a timeserie
#'
# Details:
#' see ref below
#'
# Arguments:
#'    @param timeseries a numeric vector of the observed univariate timeseries values or a numeric matrix where the first column represents the time index and the second the observed timeseries values. Use vectors/matrices with headings. If the powerspectrum is to be plotted as well, the timeseries lenght should be even number.
#'    @param winsize is the size of the rolling window expressed as percentage of the timeseries length (must be numeric between 0 and 100). Default is 50\%.
#'    @param bandwidth is the bandwidth used for the Gaussian kernel when gaussian filtering is applied. It is expressed as percentage of the timeseries length (must be numeric between 0 and 100). Alternatively it can be given by the bandwidth selector \code{\link{bw.nrd0}} (Default).
#'    @param detrending the timeseries can be detrended/filtered prior to analysis. There are four options: \code{gaussian} filtering, \code{loess} fitting, \code{linear} detrending and \code{first-differencing}. Default is \code{no} detrending.
#'    @param span parameter that controls the degree of smoothing (numeric between 0 and 100, Default 25). see more on loess{stats}
#'    @param degree the degree of polynomial to be used for when loess fitting is applied, normally 1 or 2 (Default). see more on loess{stats}
#'    @param logtransform logical. If TRUE data are logtransformed prior to analysis as log(X+1). Default is FALSE.
#'    @param interpolate logical. If TRUE linear interpolation is applied to produce a timeseries of equal length as the original. Default is FALSE (assumes there are no gaps in the timeseries).
#'    @param AR_n logical. If TRUE the best fitted AR(n) model is fitted to the data. Default is FALSE.
#'    @param powerspectrum logical. If TRUE the power spectrum within each rolling window is plotted. Default is FALSE. 
#' 
# Returns:
#'   @return \code{generic_ews} returns a matrix that contains:
#'   @return \item{tim}{the time index.}
#'   @return \item{ar1}{the \code{autoregressive coefficient ar(1)} of a first order AR model fitted on the data within the rolling window.}
#'   @return \item{sd}{the \code{standard deviation} of the data estimated within each rolling window.}
#'   @return \item{sk}{the \code{skewness} of the data estimated within each rolling window.}
#'   @return \item{kurt}{the \code{kurtosis} of the data estimated within each rolling window.}
#'   @return \item{cv}{the \code{coefficient of variation} of the data estimated within each rolling window.}
#'   @return \item{returnrate}{the return rate of the data estimated as \code{1-ar(1)} cofficient within each rolling window.}
#'   @return \item{densratio}{the \code{density ratio} of the power spectrum of the data estimated as the ratio of low frequencies over high frequencies within each rolling window.}
#'   @return \item{acf1}{the \code{autocorrelation at first lag} of the data estimated within each rolling window.}
#'
#' In addition, \code{generic_ews} returns three plots. The first plot contains the original data, the detrending/filtering applied and the residuals (if selected), and all the moment statistics. For each statistic trends are estimated by the nonparametric Kendall tau correlation.  The second plot, if asked, quantifies resilience indicators fitting AR(n) selected by the Akaike Information Criterion. The third plot, if asked, is the power spectrum estimated by \code{\link{spec.ar}} for all frequencies within each rolling window.
#'  
#' @export
#' 
#' @author Vasilis Dakos \email{vasilis.dakos@@gmail.com}
#' @references Ives, A. R. (1995). 'Measuring resilience in stochastic systems.' \emph{Ecological Monographs} 65: 217-233
#' 
#' Dakos, V., et al (2008). 'Slowing down as an early warning signal for abrupt climate change.' \emph{Proceedings of the National Academy of Sciences} 105(38): 14308-14312 
#' 
#' Dakos, V., et al (2012).'Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data.' \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010 
#' @seealso \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}}; \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}}; \code{\link{surrogates_ews}}; \code{\link{ch_ews}}; \code{\link{movpotential_ews}}; \code{\link{livpotential_ews}}; 
#'
#' @importFrom moments skewness
#' @importFrom moments kurtosis
#'
#' @examples 
#' data(foldbif)
#' out=generic_ews(foldbif,winsize=50,detrending='gaussian',
#' bandwidth=5,logtransform=FALSE,interpolate=FALSE)
#' @keywords early-warning

# Author: Vasilis Dakos, January 2, 2012

generic_ews <- function(timeseries, winsize = 50, detrending = c("no", "gaussian", 
    "loess", "linear", "first-diff"), bandwidth = NULL, span = NULL, degree = NULL, 
    logtransform = FALSE, interpolate = FALSE, AR_n = FALSE, powerspectrum = FALSE) {
    
    # timeseries<-ts(timeseries)
    timeseries <- data.matrix(timeseries)  #strict data-types the input data as tseries object for use in later steps
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
            bw <- round(length(Y) * bandwidth/100)
        }
        smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = bw, range.x = range(timeindex), 
            x.points = timeindex)
        nsmY <- Y - smYY$y
        smY <- smYY$y
    } else if (detrending == "linear") {
        nsmY <- resid(lm(Y ~ timeindex))
        smY <- fitted(lm(Y ~ timeindex))
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
    } else if (detrending == "first-diff") {
        nsmY <- diff(Y)
        timeindexdiff <- timeindex[1:(length(timeindex) - 1)]
    } else if (detrending == "no") {
        smY <- Y
        nsmY <- Y
    }
    
    
    # Rearrange data for indicator calculation
    mw <- round(length(Y) * winsize/100)
    omw <- length(nsmY) - mw + 1  ##number of moving windows
    low <- 6
    high <- omw
    nMR <- matrix(data = NA, nrow = mw, ncol = omw)
    x1 <- 1:mw
    for (i in 1:omw) {
        Ytw <- nsmY[i:(i + mw - 1)]
        nMR[, i] <- Ytw
    }
    
    # Calculate indicators
    nARR <- numeric()
    nSD <- numeric()
    nSK <- numeric()
    nKURT <- numeric()
    nACF <- numeric()
    nDENSITYRATIO <- numeric()
    nSPECT <- matrix(0, nrow = omw, ncol = ncol(nMR))
    nCV <- numeric()
    smARall <- numeric()
    smARmaxeig <- numeric()
    detB <- numeric()
    ARn <- numeric()
    
    nSD <- apply(nMR, 2, sd, na.rm = TRUE)
    for (i in 1:ncol(nMR)) {
        nYR <- ar.ols(nMR[, i], aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)
        nARR[i] <- nYR$ar
        # nSD[i]<-sapply(nMR[,i], sd, na.rm = TRUE)#sd(nMR[,i], na.rm = TRUE)
        nSK[i] <- abs(moments::skewness(nMR[, i], na.rm = TRUE))
        nKURT[i] <- moments::kurtosis(nMR[, i], na.rm = TRUE)
        nCV[i] <- nSD[i]/mean(nMR[, i])
        ACF <- acf(nMR[, i], lag.max = 1, type = c("correlation"), plot = FALSE)
        nACF[i] <- ACF$acf[2]
        spectfft <- spec.ar(nMR[, i], n.freq = omw, plot = FALSE, order = 1)
        nSPECT[, i] <- spectfft$spec
        nDENSITYRATIO[i] <- spectfft$spec[low]/spectfft$spec[high]
        
        if (AR_n) {
            ## RESILIENCE IVES 2003 Indicators based on AR(n)
            ARall <- ar.ols(nMR[, i], aic = TRUE, order.max = 6, demean = F, intercept = F)
            smARall[i] <- ARall$ar[1]
            ARn[i] <- ARall$order
            roots <- Mod(polyroot(c(rev(-ARall$ar), 1)))
            smARmaxeig[i] <- max(roots)
            detB[i] <- (prod(roots))^(2/ARn[i])
        }
    }
    
    nRETURNRATE = 1/nARR
    
    # Estimate Kendall trend statistic for indicators
    timevec <- seq(1, length(nARR))
    KtAR <- cor.test(timevec, nARR, alternative = c("two.sided"), method = c("kendall"), 
        conf.level = 0.95)
    KtACF <- cor.test(timevec, nACF, alternative = c("two.sided"), method = c("kendall"), 
        conf.level = 0.95)
    KtSD <- cor.test(timevec, nSD, alternative = c("two.sided"), method = c("kendall"), 
        conf.level = 0.95)
    KtSK <- cor.test(timevec, nSK, alternative = c("two.sided"), method = c("kendall"), 
        conf.level = 0.95)
    KtKU <- cor.test(timevec, nKURT, alternative = c("two.sided"), method = c("kendall"), 
        conf.level = 0.95)
    KtDENSITYRATIO <- cor.test(timevec, nDENSITYRATIO, alternative = c("two.sided"), 
        method = c("kendall"), conf.level = 0.95)
    KtRETURNRATE <- cor.test(timevec, nRETURNRATE, alternative = c("two.sided"), 
        method = c("kendall"), conf.level = 0.95)
    KtCV <- cor.test(timevec, nCV, alternative = c("two.sided"), method = c("kendall"), 
        conf.level = 0.95)
    
    # Plotting Generic Early-Warnings
    dev.new()
    par(mar = (c(0, 2, 0, 1) + 0), oma = c(7, 2, 3, 1), mfrow = c(5, 2))
    plot(timeindex, Y, type = "l", ylab = "", xlab = "", xaxt = "n", las = 1, xlim = c(timeindex[1], 
        timeindex[length(timeindex)]))
    if (detrending == "gaussian") {
        lines(timeindex, smY, type = "l", ylab = "", xlab = "", xaxt = "n", col = 2, 
            las = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))
    }
    if (detrending == "loess") {
        lines(timeindex, smY, type = "l", ylab = "", xlab = "", xaxt = "n", col = 2, 
            las = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))
    }
    if (detrending == "no") {
        plot(c(0, 1), c(0, 1), ylab = "", xlab = "", yaxt = "n", xaxt = "n", type = "n", 
            las = 1)
        text(0.5, 0.5, "no residuals - no detrending")
    } else if (detrending == "first-diff") {
        limit <- max(c(max(abs(nsmY))))
        plot(timeindexdiff, nsmY, ylab = "", xlab = "", type = "l", xaxt = "n", las = 1, 
            ylim = c(-limit, limit), xlim = c(timeindexdiff[1], timeindexdiff[length(timeindexdiff)]))
        legend("topleft", "first-differenced", bty = "n")
    } else {
        limit <- max(c(max(abs(nsmY))))
        plot(timeindex, nsmY, ylab = "", xlab = "", type = "h", xaxt = "n", las = 1, 
            ylim = c(-limit, limit), xlim = c(timeindex[1], timeindex[length(timeindex)]))
        legend("topleft", "residuals", bty = "n")
    }
    plot(timeindex[mw:length(nsmY)], nARR, ylab = "", xlab = "", type = "l", xaxt = "n", 
        las = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))  #3
    legend("bottomleft", paste("Kendall tau=", round(KtAR$estimate, digits = 3)), 
        bty = "n")
    legend("topleft", "ar(1)", bty = "n")
    plot(timeindex[mw:length(nsmY)], nACF, ylab = "", xlab = "", type = "l", xaxt = "n", 
        las = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))  #4
    legend("bottomleft", paste("Kendall tau=", round(KtACF$estimate, digits = 3)), 
        bty = "n")
    legend("topleft", "acf(1)", bty = "n")
    plot(timeindex[mw:length(nsmY)], nRETURNRATE, ylab = "", xlab = "", type = "l", 
        xaxt = "n", las = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))
    legend("bottomleft", paste("Kendall tau=", round(KtRETURNRATE$estimate, digits = 3)), 
        bty = "n")
    legend("topleft", "return rate", bty = "n")
    plot(timeindex[mw:length(nsmY)], nDENSITYRATIO, ylab = "", xlab = "", type = "l", 
        xaxt = "n", las = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))
    legend("bottomleft", paste("Kendall tau=", round(KtDENSITYRATIO$estimate, digits = 3)), 
        bty = "n")
    legend("topleft", "density ratio", bty = "n")
    plot(timeindex[mw:length(nsmY)], nSD, ylab = "", xlab = "", type = "l", xaxt = "n", 
        las = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))
    legend("bottomleft", paste("Kendall tau=", round(KtSD$estimate, digits = 3)), 
        bty = "n")
    legend("topleft", "standard deviation", bty = "n")
    if (detrending == "no") {
        plot(timeindex[mw:length(nsmY)], nCV, ylab = "", xlab = "", type = "l", xaxt = "n", 
            las = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))
        legend("bottomleft", paste("Kendall tau=", round(KtCV$estimate, digits = 3)), 
            bty = "n")
        legend("topleft", "coefficient of variation", bty = "n")
    } else {
        plot(0, 0, ylab = "", xlab = "", type = "n", xaxt = "n", yaxt = "n", xlim = c(0, 
            1), ylim = c(0, 1))
        text(0.5, 0.5, "no coeff var estimated - data detrended")
    }
    plot(timeindex[mw:length(nsmY)], nSK, type = "l", ylab = "", xlab = "", las = 1, 
        cex.lab = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))
    legend("topleft", "skewness", bty = "n")
    legend("bottomleft", paste("Kendall tau=", round(KtSK$estimate, digits = 3)), 
        bty = "n")
    mtext("time", side = 1, line = 2, cex = 0.8)
    plot(timeindex[mw:length(nsmY)], nKURT, type = "l", ylab = "", xlab = "", las = 1, 
        cex.lab = 1, xlim = c(timeindex[1], timeindex[length(timeindex)]))
    legend("topleft", "kurtosis", bty = "n")
    legend("bottomleft", paste("Kendall tau=", round(KtKU$estimate, digits = 3)), 
        bty = "n")
    mtext("time", side = 1, line = 2, cex = 0.8)
    mtext("Generic Early-Warnings", side = 3, line = 0.2, outer = TRUE)  #outer=TRUE print on the outer margin
    
    # Resilience Estimators based on AR(n)
    if (AR_n) {
        dev.new()
        par(mar = (c(1, 2, 0, 1) + 0.2), oma = c(4, 2, 3, 1), mfrow = c(2, 2))
        plot(timeindex[mw:length(nsmY)], ARn, type = "p", ylab = "", xlab = "", xaxt = "n", 
            cex = 0.1, las = 1, cex.axis = 0.8, xlim = c(timeindex[1], timeindex[length(timeindex)]))  #10
        legend("topleft", "AR(n)", bty = "n")
        plot(timeindex[mw:length(nsmY)], smARmaxeig, type = "l", ylab = "", xlab = "", 
            xaxt = "n", las = 1, cex.axis = 0.8, xlim = c(timeindex[1], timeindex[length(timeindex)]))
        legend("topleft", "max eigenvalue", bty = "n")
        plot(timeindex[mw:length(nsmY)], detB, type = "l", ylab = "", xlab = "", 
            cex.lab = 1, las = 1, cex.axis = 0.8, xlim = c(timeindex[1], timeindex[length(timeindex)]))
        mtext("time", side = 1, line = 2, cex = 0.8)
        legend("topleft", "geometric mean root AR(n)", bty = "n")
        plot(timeindex[mw:length(nsmY)], smARall, type = "l", ylab = "", xlab = "", 
            cex.lab = 1, las = 1, cex.axis = 0.8, xlim = c(timeindex[1], timeindex[length(timeindex)]))
        mtext("time", side = 1, line = 2, cex = 0.8)
        legend("topleft", "b1 of AR(n)", bty = "n")
        mtext("Resilience Estimators based on AR(n)", side = 3, line = 0.2, outer = TRUE)
    }
    
    # Power spectrum
    if (powerspectrum) {
        dev.new()
        par(mar = (c(4.6, 4, 0.5, 2) + 0.2), oma = c(0.5, 1, 2, 1))
        image(x = (spectfft$freq[2:length(spectfft$freq)]), y = (seq(1, ncol(nSPECT), 
            by = 1)), log(nSPECT[2:length(spectfft$freq), ]), ylab = "rolling window", 
            xlab = "frequency", log = "x", xlim = c(spectfft$freq[2], spectfft$freq[length(spectfft$freq)]), 
            col = topo.colors(20), xaxs = "i")
        contour(x = (spectfft$freq[2:length(spectfft$freq)]), y = (seq(1, ncol(nSPECT), 
            by = 1)), log(nSPECT[2:length(spectfft$freq), ]), add = TRUE)
        mtext("Power spectrum within rolling windows", side = 3, line = 0.2, outer = TRUE)
    }
    
    # Output
    out <- data.frame(timeindex[mw:length(nsmY)], nARR, nSD, nSK, nKURT, nCV, nRETURNRATE, 
        nDENSITYRATIO, nACF)
    colnames(out) <- c("timeindex", "ar1", "sd", "sk", "kurt", "cv", "returnrate", 
        "densratio", "acf1")
    return(out)
} 
