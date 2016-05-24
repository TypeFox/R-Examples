# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Prediction loss
#' 
#' Compute the prediction loss of a model.  
#' 
#' \code{mspe} and \code{rmspe} compute the mean squared prediction error and 
#' the root mean squared prediction error, respectively.  In addition, 
#' \code{mape} returns the mean absolute prediction error, which is somewhat 
#' more robust.
#' 
#' Robust prediction loss based on trimming is implemented in \code{tmspe} and 
#' \code{rtmspe}.  To be more precise, \code{tmspe} computes the trimmed mean 
#' squared prediction error and \code{rtmspe} computes the root trimmed mean 
#' squared prediction error.  A proportion of the largest squared differences 
#' of the observed and fitted values are thereby trimmed.
#' 
#' Standard errors can be requested via the \code{includeSE} argument.  Note that 
#' standard errors for \code{tmspe} are based on a winsorized standard 
#' deviation.  Furthermore, standard errors for \code{rmspe} and \code{rtmspe} 
#' are computed from the respective standard errors of \code{mspe} and 
#' \code{tmspe} via the delta method.
#' 
#' @rdname cost
#' @name cost
#' 
#' @param y  a numeric vector or matrix giving the observed values.
#' @param yHat  a numeric vector or matrix of the same dimensions as \code{y} 
#' giving the fitted values.
#' @param trim  a numeric value giving the trimming proportion (the default is 
#' 0.25).
#' @param includeSE  a logical indicating whether standard errors should be computed 
#' as well.
#' 
#' @return If standard errors are not requested, a numeric value giving the 
#' prediction loss is returned.
#' 
#' Otherwise a list is returned, with the first component containing the 
#' prediction loss and the second component the corresponding standard error.
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Tukey, J.W. and McLaughlin, D.H. (1963) Less vulnerable confidence and 
#' significance procedures for location based on a single sample: 
#' Trimming/winsorization.  \emph{Sankhya: The Indian Journal of Statistics, 
#' Series A}, \bold{25}(3), 331--352
#' 
#' Oehlert, G.W. (1992) A note on the delta method.  \emph{The American 
#' Statistician}, \bold{46}(1), 27--29.
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvTuning}}
#' 
#' @examples 
#' # fit an MM-regression model
#' data("coleman")
#' fit <- lmrob(Y~., data=coleman)
#' 
#' # compute the prediction loss from the fitted values
#' # (hence the prediction loss is underestimated in this simple 
#' # example since all observations are used to fit the model)
#' mspe(coleman$Y, predict(fit))
#' rmspe(coleman$Y, predict(fit))
#' mape(coleman$Y, predict(fit))
#' tmspe(coleman$Y, predict(fit), trim = 0.1)
#' rtmspe(coleman$Y, predict(fit), trim = 0.1)
#' 
#' # include standard error
#' mspe(coleman$Y, predict(fit), includeSE = TRUE)
#' rmspe(coleman$Y, predict(fit), includeSE = TRUE)
#' mape(coleman$Y, predict(fit), includeSE = TRUE)
#' tmspe(coleman$Y, predict(fit), trim = 0.1, includeSE = TRUE)
#' rtmspe(coleman$Y, predict(fit), trim = 0.1, includeSE = TRUE)
#' 
#' @keywords utilities

NULL

## mean squared prediction error
#' @rdname cost
#' @export
mspe <- function(y, yHat, includeSE = FALSE) {
    residuals2 <- (y - yHat)^2  # squared residuals
    if(!is.null(dim(y))) {
        residuals2 <- rowSums(residuals2)  # squared norm in multivariate case
    }
    res <- mean(residuals2)
    if(isTRUE(includeSE)) {
        res <- list(mspe=res, se=sd(residuals2)/sqrt(nobs(residuals2)))
    }
    res
}

## root mean squared prediction error
#' @rdname cost
#' @export
rmspe <- function(y, yHat, includeSE = FALSE) {
    includeSE <- isTRUE(includeSE)
    res <- mspe(y, yHat, includeSE=includeSE)
    if(includeSE) {
        rmspe <- sqrt(res$mspe)
        res <- list(rmspe=rmspe, se=res$se/(2*rmspe))
    } else res <- sqrt(res)
    res
}

## mean absolute prediction error
#' @rdname cost
#' @export
mape <- function(y, yHat, includeSE = FALSE) {
    absResiduals <- abs(y - yHat)  # absolue residuals
    if(!is.null(dim(y))) {
        absResiduals <- rowSums(absResiduals)  # norm in multivariate case
    }
    res <- mean(absResiduals)
    if(isTRUE(includeSE)) {
        res <- list(mape=res, se=sd(absResiduals)/sqrt(nobs(absResiduals)))
    }
    res
}

## trimmed mean squared prediction error
#' @rdname cost
#' @export
tmspe <- function(y, yHat, trim = 0.25, includeSE = FALSE) {
    n <- nobs(y)
    h <- n - floor(n*trim)
    residuals2 <- (y - yHat)^2  # squared residuals
    if(!is.null(dim(y))) {
        residuals2 <- rowSums(residuals2)  # squared norm in multivariate case
    }
    residuals2 <- sort(residuals2)       # sort squared residuals
    res <- mean(residuals2[seq_len(h)])  # mean over h smallest values
    if(isTRUE(includeSE)) {
        # standard error of the trimmed mean is based on a winsorized 
        # standard deviation
        alpha <- 1 - trim
        if(h < n) {
            q <- quantile(residuals2, alpha)  # quantile for winsorization
            residuals2[(h+1):n] <- q          # replace tail with quantile
        }
        res <- list(tmspe=res, se=sd(residuals2)/(alpha*sqrt(n)))
    }
    res
}

## root trimmed mean squared prediction error
#' @rdname cost
#' @export
rtmspe <- function(y, yHat, trim = 0.25, includeSE = FALSE) {
    includeSE <- isTRUE(includeSE)
    res <- tmspe(y, yHat, trim=trim, includeSE=includeSE)
    if(includeSE) {
        rtmspe <- sqrt(res$tmspe)
        res <- list(rtmspe=rtmspe, se=res$se/(2*rtmspe))
    } else res <- sqrt(res)
    res
}
