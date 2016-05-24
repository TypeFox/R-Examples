# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

#' Hill estimator
#' 
#' The Hill estimator uses the maximum likelihood principle to estimate the
#' shape parameter of a Pareto distribution.
#' 
#' The arguments \code{k} and \code{x0} of course correspond with each other.
#' If \code{k} is supplied, the threshold \code{x0} is estimated with the \eqn{n
#' - k} largest value in \code{x}, where \eqn{n} is the number of observations.
#' On the other hand, if the threshold \code{x0} is supplied, \code{k} is given
#' by the number of observations in \code{x} larger than \code{x0}.  Therefore,
#' either \code{k} or \code{x0} needs to be supplied.  If both are supplied,
#' only \code{k} is used (mainly for back compatibility).
#' 
#' @param x a numeric vector.
#' @param k the number of observations in the upper tail to which the Pareto
#' distribution is fitted.
#' @param x0 the threshold (scale parameter) above which the Pareto distribution
#' is fitted.
#' @param w an optional numeric vector giving sample weights.
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
#' \code{\link{thetaPDC}}, \code{\link{thetaWML}}, \code{\link{thetaISE}},
#' \code{\link{minAMSE}}
#' 
#' @references Hill, B.M. (1975) A simple general approach to inference about
#' the tail of a distribution. \emph{The Annals of Statistics}, \bold{3}(5),
#' 1163--1174.
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
#' thetaHill(eusilc$eqIncome, k = ts$k, w = eusilc$db090)
#' 
#' # using threshold
#' thetaHill(eusilc$eqIncome, x0 = ts$x0, w = eusilc$db090)
#' 
#' @export

thetaHill <- function (x, k = NULL, x0 = NULL, w = NULL) {
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
    .thetaHill(x, k, x0, w)
}

.thetaHill <- function (x, k = NULL, x0 = NULL, w = NULL) {
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
    if(haveW) {
        ## return weighted Hill estimate
        w <- w[(n-k+1):n]
        sum(w)/sum(w*(log(x[(n-k+1):n]) - log(x0)))
    } else {
        ## return Hill estimate
#        k/(sum(log(x[(n-k+1):n])) - k*log(x0))
        k/sum(log(x[(n-k+1):n]) - log(x0))  # should be numerically more stable
    }
}
