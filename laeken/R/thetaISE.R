# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

#' Integrated squared error (ISE) estimator
#' 
#' The integrated squared error (ISE) estimator estimates the shape parameter of
#' a Pareto distribution based on the relative excesses of observations above a
#' certain threshold.
#' 
#' The arguments \code{k} and \code{x0} of course correspond with each other.
#' If \code{k} is supplied, the threshold \code{x0} is estimated with the \eqn{n
#' - k} largest value in \code{x}, where \eqn{n} is the number of observations.
#' On the other hand, if the threshold \code{x0} is supplied, \code{k} is given
#' by the number of observations in \code{x} larger than \code{x0}.  Therefore,
#' either \code{k} or \code{x0} needs to be supplied.  If both are supplied,
#' only \code{k} is used (mainly for back compatibility).
#' 
#' The ISE estimator minimizes the integrated squared error (ISE) criterion with
#' a complete density model.  The minimization is carried out using %
#' \code{\link[stats]{nlm}}.  By default, the starting value is obtained % with
#' the Hill estimator (see \code{\link{thetaHill}}).
#' \code{\link[stats]{optimize}}.
#' 
#' @param x a numeric vector.
#' @param k the number of observations in the upper tail to which the Pareto
#' distribution is fitted.
#' @param x0 the threshold (scale parameter) above which the Pareto distribution
#' is fitted.
#' @param w an optional numeric vector giving sample weights.
#' @param \dots additional arguments to be passed to
#' \code{\link[stats]{optimize}} (see \dQuote{Details}).
#' 
#' @return The estimated shape parameter.
#' 
#' @note The arguments \code{x0} for the threshold (scale parameter) of the
#' Pareto distribution and \code{w} for sample weights were introduced in
#' version 0.2.
#' 
#' @author Andreas Alfons and Josef Holzer
#' 
#' @seealso \code{\link{paretoTail}}, \code{\link{fitPareto}},
#' \code{\link{thetaPDC}}, \code{\link{thetaHill}}
#' 
#' @references 
#' A. Alfons and M. Templ (2013) Estimation of Social Exclusion Indicators 
#' from Complex Surveys: The \R Package \pkg{laeken}.  \emph{Journal of 
#' Statistical Software}, \bold{54}(15), 1--25.  URL 
#' \url{http://www.jstatsoft.org/v54/i15/}
#' 
#' A. Alfons, M. Templ, P. Filzmoser (2013) Robust estimation of economic 
#' indicators from survey samples based on Pareto tail modeling. \emph{Journal 
#' of the Royal Statistical Society, Series C}, \bold{62}(2), 271--286.
#' 
#' Vandewalle, B., Beirlant, J., Christmann, A., and Hubert, M.
#' (2007) A robust estimator for the tail index of Pareto-type 
#' distributions.  \emph{Computational Statistics & Data Analysis}, 
#' \bold{51}(12), 6252--6268.
#' 
#' @keywords manip
#' 
#' @examples
#' data(eusilc)
#' # equivalized disposable income is equal for each household
#' # member, therefore only one household member is taken
#' eusilc <- eusilc[!duplicated(eusilc$db030),]
#' 
#' # estimate threshold
#' ts <- paretoScale(eusilc$eqIncome, w = eusilc$db090)
#' 
#' # using number of observations in tail
#' thetaISE(eusilc$eqIncome, k = ts$k, w = eusilc$db090)
#' 
#' # using threshold
#' thetaISE(eusilc$eqIncome, x0 = ts$x0, w = eusilc$db090)
#' 
#' @export

thetaISE <- function(x, k = NULL, x0 = NULL, w = NULL, ...) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    haveK <- !is.null(k)
    if(haveK) {  # if 'k' is supplied, it is always used
        if(!is.numeric(k) || length(k) == 0 || k[1] < 1) {
            stop("'k' must be a positive integer")
        } else k <- k[1]
    } else if(!is.null(x0)) {  # otherwise 'x0' (threshold) is used
        if(!is.numeric(x0) || length(x0) == 0) stop("'x0' must be numeric")
        else x0 <- x0[1]
    } else stop("either 'k' or 'x0' must be supplied")
    haveW <- !is.null(w)
    if(haveW) {  # sample weights are supplied
        if(!is.numeric(w) || length(w) != length(x)) {
            stop("'w' must be numeric vector of the same length as 'x'")
        }
        if(any(w < 0)) stop("negative weights in 'w'")
        if(any(i <- is.na(x))) {  # remove missing values
            x <- x[!i]
            w <- w[!i]
        }
        # sort values and sample weights
        order <- order(x)
        x <- x[order]
        w <- w[order]
    } else {  # no sample weights
        if(any(i <- is.na(x))) x <- x[!i]  # remove missing values
        x <- sort(x)  # sort values
    }
    .thetaISE(x, k, x0, w, ...)
}

# internal function that assumes that data are ok and sorted
.thetaISE <- function(x, k = NULL, x0 = NULL, w = NULL, 
        tol = .Machine$double.eps^0.25, ...) {
    n <- length(x)  # number of observations
    haveK <- !is.null(k)
    haveW <- !is.null(w)
    if(haveK) {  # 'k' is supplied, threshold is determined
        if(k >= n) stop("'k' must be smaller than the number of observed values")
        x0 <- x[n-k]  # threshold (scale parameter)
    } else {  # 'k' is not supplied, it is determined using threshold
        # values are already sorted
        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
        k <- length(which(x > x0))
    }
    ## computations
    y <- x[(n-k+1):n]/x0  # relative excesses
    if(haveW) {
        wTail <- w[(n-k+1):n]
        ## weighted integrated squared error distance criterion
        # w ... sample weights
        ISE <- function(theta, y, w) {
            f <- theta*y^(-1-theta)
            weighted.mean(theta^2/(2*theta+1) - 2*f, w)
        }
    } else {
        wTail <- NULL
        ## integrated squared error distance criterion
        # w ... sample weights (not needed here, only available to have the 
        #       same function definition)
        ISE <- function(theta, y, w) {
            f <- theta*y^(-1-theta)
            mean(theta^2/(2*theta+1) - 2*f)
        }
    }
    ## optimize
    localOptimize <- function(f, interval = NULL, tol, ...) {
        if(is.null(interval)) {
            p <- if(haveK) .thetaHill(x, k, w=w) else .thetaHill(x, x0=x0, w=w)
            interval <- c(0 + tol, 3 * p)  # default interval
        }
        optimize(f, interval, ...)
    }
    localOptimize(ISE, y=y, w=wTail, tol=tol, ...)$minimum
}
