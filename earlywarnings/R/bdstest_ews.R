# BDS test Early Warning Signals Author: Stephen R Carpenter, 22 Oct 2011
# Modified by: Vasilis Dakos, January 1, 2012

#' @importFrom tseries bds.test

BDSboot <- function(X, varname, nboot, epsvec, emb) {
    # begin function
    
    StdEpsAll <- X  # name of variable for BDS
    neps <- length(epsvec)
    # Compute and print BDS test
    print("***********************************************", quote = FALSE)
    print(c("BDS test for ", varname), quote = FALSE)
    print(c("Embedding dimension = ", emb), quote = FALSE)
    BDS.data <- bds.test(StdEpsAll, m = emb, epsvec)
    print("BDS statistics for Nominal Data at each Epsilon", quote = FALSE)
    print(round(BDS.data$statistic, 3))
    print("P value based on standard normal", quote = FALSE)
    print(round(BDS.data$p.value, 3))
    # Bootstrap the BDS test
    nobs <- length(StdEpsAll)
    bootmat <- matrix(0, nrow = emb - 1, ncol = neps)  # matrix to count extreme BDS values
    for (i in 1:nboot) {
        # start bootstrap loop
        epsboot <- sample(StdEpsAll, nobs, replace = TRUE)
        BDS.boot <- bds.test(epsboot, m = emb, epsvec)
        for (im in 1:(emb - 1)) {
            # loop over embedding dimensions
            bootvec <- BDS.boot$statistic[im, ]
            N.above <- ifelse(bootvec > BDS.data$statistic[im, ], 1, 0)
            bootmat[im, ] <- bootmat[im, ] + N.above
        }
        # Report progress: if hash is removed from the next two lines, the program will
        # report each time an iteration is completed cat('iteration = ',i,' of
        # ',nboot,'\n') # flush.console()
    }  # end bootstrap loop
    print(" ", quote = FALSE)
    print(c("Bootstrap P estimates for ", varname), quote = FALSE)
    print(c("Bootstrap iterations = ", nboot), quote = FALSE)
    Pboot <- bootmat/nboot
    for (im in 1:(emb - 1)) {
        print(c("For embedding dimension =", im + 1), quote = FALSE)
        print(c("For epsilon = ", round(epsvec, 3), "bootstrap P = "), quote = FALSE)
        print(Pboot[im, ])
    }
    print("**********************************************************", quote = FALSE)
}  # end function


#' Description: BDS test Early Warning Signals
#'
#' \code{bdstest_ews} is used to estimate the BDS statistic to detect nonlinearity in the residuals of a timeseries after first-difference detrending, fitting an ARMA(p,q) model, and fitting a GARCH(0,1) model. The function is making use of \code{bds.test} from the tseries package.
#'
# Details:
#' See also \code{bds.test{tseries}} for more details. The function requires the installation of packages \code{tseries} and \code{quadprog} that are not available under Linux and need to be manually installed under Windows.
#'
#' Arguments:
#'    @param timeseries a numeric vector of the observed univariate timeseries values or a numeric matrix where the first column represents the time index and the second the observed timeseries values. Use vectors/matrices with headings.
#'    @param ARMAoptim is the order of the \code{ARMA(p,q)} model to be fitted on the original timeseries. If TRUE the best ARMA model based on AIC is applied. If FALSE the \code{ARMAorder} is used.
#'    @param ARMAorder is the order of the \code{AR(p)} and \code{MA(q)} process to be fitted on the original timeseries. Default is \code{p=1} \code{q=0}.
#'    @param GARCHorder fits a GARCH model on the original timeseries where \code{GARCHorder[1]} is the GARCH part and \code{GARCHorder[2]} is the ARCH part.
#'    @param embdim is the embedding dimension (2, 3,... \code{embdim}) up to which the BDS test will be estimated (must be numeric). Default value is 3.
#'    @param epsilon is a numeric vector that is used to scale the standard deviation of the timeseries. The BDS test is computed for each element of epsilon. Default is 0.5, 0.75 and 1.
#'    @param boots is the number of bootstraps performed to estimate significance p values for the BDS test. Default is 1000.
#'    @param logtransform logical. If TRUE data are logtransformed prior to analysis as log(X+1). Default is FALSE.
#'    @param interpolate logical. If TRUE linear interpolation is applied to produce a timeseries of equal length as the original. Default is FALSE (assumes there are no gaps in the timeseries).
#' 
# Values:
#' @return \code{bdstest_ews} returns output on the R console that summarizes the BDS test statistic for all embedding dimensions and  \code{epsilon} values used, and for first-differenced data, ARMA(p.q) residuals, and GARCH(0,1) residuals). Also the significance p values are returned estimated both by comparing to a standard normal distribution and by bootstrapping.
#' 
#' In addition, \code{bdstest_ews} returns a plot with the original timeseries, the residuals after first-differencing, and fitting the ARMA(p,q) and GARCH(0,1) models. Also the autocorrelation \code{\link{acf}} and partial autocorrelation \code{\link{pacf}} functions are estimated serving as guides for the choice of lags of the linear models fitted to the data.
#'  
#' @export
#' 
#' @author S. R. Carpenter, modified by V. Dakos
#' @references J. B. Cromwell, W. C. Labys and M. Terraza (1994): Univariate Tests for Time Series Models, Sage, Thousand Oaks, CA, pages 32-36.
#' 
#' Dakos, V., et al (2012).'Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data.' \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010 
#' @seealso 
#' \code{\link{generic_ews}}; 
#' \code{\link{ddjnonparam_ews}}; 
#' \code{\link{bdstest_ews}}; 
#' \code{\link{sensitivity_ews}}; 
#' \code{\link{surrogates_ews}}; 
#' \code{\link{ch_ews}}; 
#' \code{\link{movpotential_ews}}; 
#' \code{\link{livpotential_ews}}
# ; \code{\link{timeVAR_ews}}; \code{\link{thresholdAR_ews}}
#' @importFrom tseries garch
#' @import quadprog
#' @examples #
#' data(foldbif)
#' bdstest_ews(foldbif,ARMAoptim=FALSE,ARMAorder=c(1,0),embdim=3,epsilon=0.5, boots=200,logtransform=FALSE,interpolate=FALSE)
#' @keywords early-warning
#' 
bdstest_ews <- function(timeseries, ARMAoptim = TRUE, ARMAorder = c(1, 0), GARCHorder = c(0, 
    1), embdim = 3, epsilon = c(0.5, 0.75, 1), boots = 1000, logtransform = FALSE, 
    interpolate = FALSE) {
    
    timeseries <- ts(timeseries)  #strict data-types the input data as tseries object for use in later steps
    if (ncol(timeseries) == 1) {
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
    
    # Detrend the data
    Eps1 <- diff(Y)
    
    # Define BDS parameters
    nboot <- boots
    emb <- embdim  # embedding dimension
    eps.sd <- sd(as.vector(Eps1))
    epsvec <- round(eps.sd * epsilon, 6)
    
    # Run BDS with bootstrapping
    BDSboot(Eps1, c("Detrended data"), nboot, epsvec, emb)
    
    # Fit ARMA model based on AIC
    if (ARMAoptim == TRUE) {
        arma = matrix(, 4, 5)
        for (ij in 1:4) {
            for (jj in 0:4) {
                ARMA <- arima(Y, order = c(ij, 0, jj), method = "ML", include.mean = TRUE)
                arma[ij, jj + 1] = ARMA$aic
                ind = which(arma == min(arma), arr.ind = TRUE)
                armafit <- arima(Y, order = c(ind[1], 0, ind[2] - 1), include.mean = TRUE)
                print("ARMA model", quote = FALSE)
                print(armafit, digits = 4)
                Eps2 <- armafit$residuals  #[2:armafit$n.used]
            }
        }
    } else {
        armafit <- arima(Y, order = c(ARMAorder[1], 0, ARMAorder[2]), include.mean = TRUE)
        print("ARMA model", quote = FALSE)
        print(armafit, digits = 4)
        Eps2 <- armafit$residuals  #[2:armafit$n.used]
    }
    
    # Define BDS parameters
    nboot <- boots
    emb <- embdim  # embedding dimension
    eps.sd <- sd(as.vector(Eps2))
    epsvec <- round(eps.sd * epsilon, 6)
    
    # Run BDS with bootstrapping
    BDSboot(Eps2, c("ARMA model residuals"), nboot, epsvec, emb)
    
    # Fit GARCH(0,1) model to detrended data
    Gfit <- garch(Y, order = c(GARCHorder[1], GARCHorder[2]))
    print("GARCH(0,1) model fit to detrended data", quote = FALSE)
    print(Gfit, digits = 4)
    
    Eps3 <- Gfit$residuals[2:length(Y)]
    
    # Define BDS parameters
    nboot <- boots
    emb <- embdim  # embedding dimension
    eps.sd <- sd(as.vector(Eps3))
    epsvec <- round(eps.sd * epsilon, 6)
    
    # Run BDS with bootstrapping
    BDSboot(Eps3, c("GARCH(0,1) model residuals"), nboot, epsvec, emb)
    
    # Plot the data
    dev.new()
    par(fig = c(0, 0.5, 0.5, 1), mar = c(3, 4, 3, 2), cex.axis = 0.8, cex.lab = 1, 
        mfrow = c(2, 2), mgp = c(1.5, 0.5, 0), oma = c(1, 1, 2, 1))
    plot(timeindex, Y, type = "l", col = "black", lwd = 1.7, xlab = "time", ylab = "data")
    par(fig = c(0.5, 1, 0.8, 0.95), mar = c(0, 4, 0, 2), new = TRUE)
    plot(timeindex[1:(length(timeindex) - 1)], Eps1, type = "l", col = "red", lwd = 1, 
        xlab = "", ylab = "", xaxt = "n")
    legend("topright", "first-diff", bty = "n", cex = 0.8)
    par(fig = c(0.5, 1, 0.65, 0.8), new = TRUE)
    plot(timeindex[1:(length(timeindex))], Eps2, type = "l", xlab = "", ylab = "residuals", 
        xaxt = "n", col = "green", lwd = 1)
    legend("topright", "AR", bty = "n", cex = 0.8)
    par(fig = c(0.5, 1, 0.5, 0.65), new = TRUE)
    plot(timeindex[1:(length(timeindex) - 1)], Eps3, type = "l", col = "blue", xlab = "time", 
        ylab = "", lwd = 1)
    legend("topright", "GARCH", bty = "n", cex = 0.8)
    mtext("time", side = 1, outer = FALSE, line = 1.4, cex = 0.8)
    # Time series diagnostics
    par(fig = c(0, 0.5, 0, 0.4), mar = c(3, 4, 0, 2), new = TRUE)
    acf(Y, lag.max = 25, main = "")
    par(fig = c(0.5, 1, 0, 0.4), new = TRUE)
    pacf(Y, lag.max = 25, main = "")
    mtext("BDS_test Diagnostics", side = 3, line = 0.2, outer = TRUE)
} 
