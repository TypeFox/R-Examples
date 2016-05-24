#' Unit root test for events considering a structrual break
#' 
#' This function performs the Zivot & Andrews unit root test, which allows a break at an unknown point in either the intercept, 
#' the linear trend or in both.

#' 
#' @details This function is written refering to the \code{ur.za} function
#' in the \pkg{urza} package (Pfaff 2008), but it speeds up executation using the \pkg{RcppArmadillo} package. Allowing a structrual break, 
#' this function returns flag to be 0 if the time series  is stationary and 1 if it is a unit root process.

#' 
#' @param y a vector or a time series.
#' @param model Three choices: ``intercept', ``trend' or ``both'.
#' @param lag a scalar chosen as lag.
#' @return a list consisting of:
#' \item{flag}{0 if the time series is 
#' is stationary; 1 if it is a unit root process.}
#' \item{teststat}{ZA unit root test statistic.} 
#' @seealso \code{\link{noiseTests}}
#' 
#' @references Eric Zivot and Donald W K Andrews (1992). Further evidence on the great crash, the oil-price shock, and the unit-root hypothesis. \emph{Journal of Business & Economic Statistics}, 
#' \bold{20}(1), 25-44. \url{http://dx.doi.org/10.1198/073500102753410372}.
#' 
#' @references Pfaff, Bernhard (2008). Analysis of Integrated and Cointegrated Time Series with R. Second Edition. Springer, New York.
#' \url{http://www.springer.com/statistics/statistical+theory+and+methods/book/978-0-387-75966-1}.

#' @export
#' @examples
#' # this is a box function
#' set.seed(123)
#' x=cbfs_red('box')
#' ur.za.fast(x,'both')
#' # this is a cliff-ramp
#' set.seed(123)
#' x=cbfs_red('cr')
#' ur.za.fast(x,'both')
#' # this is a random walk process
#' set.seed(123)
#' x=cumsum(rnorm(300))
#' ur.za.fast(x,'both')


ur.za.fast <- function(y, model = c("intercept", "trend", "both"), lag = NULL) {
    n <- length(y)
    model <- match.arg(model)
    if (is.null(lag)) 
        lag <- 0
    lag <- as.integer(lag)
    datmat <- matrix(NA, n, lag + 4)
    idx <- 2:(n - 2)
    trend <- seq(1, n)
    datmat[, 1:4] <- cbind(y, 1, c(NA, y)[1:n], trend)
    colnames(datmat)[1:4] <- c("y", "intercept", "y.l1", "trend")
    if (lag > 0) {
        for (i in 1:lag) {
            datmat[, i + 3] <- c(rep(NA, i + 1), diff(y))[1:n]
        }
        colnames(datmat) <- c("y", "y.l1", "trend", paste("y.dl", 1:lag, sep = ""))
    }
    if (model == "intercept") {
        roll <- function(z) {
            du <- c(rep(0, z), rep(1, (n - z)))
            rollmat <- cbind(datmat, du)
            roll.reg <- tryCatch(fastLmPure(rollmat[2:dim(rollmat)[1], 2:dim(rollmat)[2]], rollmat[2:dim(rollmat)[1], 1]),error=function(e){1})
            if ( class(roll.reg)!='list') {
                roll.reg <- coef(summary(lm(as.data.frame(rollmat))))
                (roll.reg[2, 1] - 1)/roll.reg[2, 2]
            } else {
                (roll.reg$coefficients[2] - 1)/roll.reg$stderr[2]
            }
        }
        roll.stat <- sapply(idx, roll)
        cval <- c(-5.34, -4.8, -4.58)
        bpoint <- which.min(roll.stat)
    } else if (model == "trend") {
        roll <- function(z) {
            dt <- c(rep(0, z), 1:(n - z))
            rollmat <- cbind(datmat, dt)
            roll.reg <- tryCatch(fastLmPure(rollmat[2:dim(rollmat)[1], 2:dim(rollmat)[2]], 
                                            rollmat[2:dim(rollmat)[1], 1]),error=function(e){1})
            if ( class(roll.reg)!='list') {
                roll.reg <- coef(summary(lm(as.data.frame(rollmat))))
                (roll.reg[2, 1] - 1)/roll.reg[2, 2]
            } else {
                (roll.reg$coefficients[2] - 1)/roll.reg$stderr[2]
            }
        }
        roll.stat <- sapply(idx, roll)
        cval <- c(-4.93, -4.42, -4.11)
        bpoint <- which.min(roll.stat)
    } else if (model == "both") {
        roll <- function(z) {
            du <- c(rep(0, z), rep(1, (n - z)))
            dt <- c(rep(0, z), 1:(n - z))
            rollmat <- cbind(datmat, du, dt)
            roll.reg <- tryCatch(fastLmPure(rollmat[2:dim(rollmat)[1], 2:dim(rollmat)[2]], rollmat[2:dim(rollmat)[1], 1]),error=function(e){1})
            if ( class(roll.reg)!='list') {
                roll.reg <- coef(summary(lm(as.data.frame(rollmat))))
                (roll.reg[2, 1] - 1)/roll.reg[2, 2]
            } else {
                (roll.reg$coefficients[2] - 1)/roll.reg$stderr[2]
            }
        }
        roll.stat <- sapply(idx, roll)
        cval <- c(-5.57, -5.08, -4.82)
        bpoint <- which.min(roll.stat)
    }
    teststat <- roll.stat[bpoint]
    flag = ifelse(teststat < cval[2], 0, 1)
    # results = list(teststat = teststat, cval = cval, bpoint = bpoint,flag=flag)
    return(list(flag = flag, teststat = teststat))
} 
