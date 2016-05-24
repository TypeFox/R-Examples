# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

#' Estimate the scale parameter of a Pareto distribution
#'
#' Estimate the scale parameter of a Pareto distribution, i.e., the threshold
#' for Pareto tail modeling.
#'
#' Van Kerm's formula is given by \deqn{\min(\max(2.5 \bar{x}, q(0.98),
#' q(0.97))),}{min(max(2.5 m(x), q(0.98)), q(0.97)),} where \eqn{\bar{x}}{m(x)}
#' denotes the weighted mean and \eqn{q(.)} denotes weighted quantiles.  This
#' function allows to compute generalizations of Van Kerm's formula, where the
#' mean can be replaced by the median and different quantiles can be used.
#'
#' @aliases print.paretoScale
#'
#' @param x a numeric vector.
#' @param w an optional numeric vector giving sample weights.
#' @param groups an optional vector or factor specifying groups of elements of
#' \code{x} (e.g., households).  If supplied, each group of observations is
#' expected to have the same value in \code{x} (e.g., household income).  Only
#' the values of every first group member to appear are used for estimating the
#' threshold (scale parameter).
#' @param method a character string specifying the estimation method.  If
#' \code{"VanKerm"}, Van Kerm's method is used, which is a rule of thumb
#' specifically designed for the equivalized disposable income in EU-SILC data
#' (currently the only method implemented).
#' @param center a character string specifying the estimation method for the
#' center of the distribution.  Possible values are \code{"mean"} for the
#' weighted mean and \code{"median"} for the weighted median.  This is used if
#' \code{method} is \code{"VanKerm"} (currently the only method implemented).
#' @param probs a numeric vector of length two giving probabilities to be used
#' for computing weighted quantiles of the distribution.  Values should be close
#' to 1 such that the quantiles correspond to the upper tail.  This is used if
#' \code{method} is \code{"VanKerm"} (currently the only method implemented).
#' @param na.rm a logical indicating whether missing values in \code{x} should
#' be omitted.
#'
#' @returnClass paretoScale
#' @returnItem x0 the threshold (scale parameter).
#' @returnItem k the number of observations in the tail (i.e., larger than the
#' threshold).
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{minAMSE}}, \code{\link{paretoQPlot}},
#' \code{\link{meanExcessPlot}}
#'
#' @references
#' A. Alfons and M. Templ (2013) Estimation of Social Exclusion Indicators
#' from Complex Surveys: The \R Package \pkg{laeken}.  \emph{Journal of
#' Statistical Software}, \bold{54}(15), 1--25.  URL
#' \url{http://www.jstatsoft.org/v54/i15/}
#'
#' Van Kerm, P. (2007) Extreme incomes and the estimation of poverty and
#' inequality indicators from EU-SILC. IRISS Working Paper Series 2007-01,
#' CEPS/INSTEAD.
#'
#' @keywords manip
#'
#' @examples
#' data(eusilc)
#' paretoScale(eusilc$eqIncome, eusilc$db090, groups = eusilc$db030)
#'
#' @export

paretoScale <- function(x, w = NULL, groups = NULL,
        method = "VanKerm", center = c("mean", "median"),
		probs = c(0.97, 0.98), na.rm = FALSE) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    useW <- !is.null(w)
    if(useW && (!is.numeric(w) || length(w) != length(x))) {
        stop("'w' must be numeric vector of the same length as 'x'")
    }
    haveGroups <- !is.null(groups)
    if(haveGroups) {
        if(!is.vector(groups) && !is.factor(groups)) {
            stop("'groups' must be a vector or factor")
        }
        if(length(groups) != length(x)) {
            stop("'groups' must have the same length as 'x'")
        }
        if(any(is.na(groups))) stop("'groups' contains missing values")
        unique <- !duplicated(groups)
        x <- x[unique]
        if(useW) w <- w[unique]
    }
#	method <- match.arg(method)  # only van Kerm's method currently implemented
	center <- match.arg(center)
	probs <- rep(probs, length.out=2)
    na.rm <- isTRUE(na.rm)
    # estimate threshold with van Kerm's method
    if(center == "mean") {
		mu <- weightedMean(x, w, na.rm=na.rm)
		q <- weightedQuantile(x, w, probs=probs, na.rm=na.rm)
	} else {
		q <- weightedQuantile(x, w, probs=c(0.5, probs), na.rm=na.rm)
		mu <- q[1]
		q <- q[-1]
	}
    x0 <- max(min(2.5*mu, q[2]), q[1])
    res <- list(x0=x0, k=length(which(x > x0)))
    class(res) <- "paretoScale"
    res
}


## print method for class "paretoScale"
#' @export
print.paretoScale <- function(x, ...) {
    cat("Threshold: ")
    cat(x$x0, ...)
    cat("\nNumber of observations in the tail: ")
    cat(x$k, ...)
    cat("\n")
}
