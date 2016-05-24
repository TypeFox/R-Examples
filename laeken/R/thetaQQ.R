# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

#' QQ-estimator
#' 
#' Estimate the shape parameter of a Pareto distribution using a
#' quantile-quantile approach.
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
#' 
#' @return The estimated shape parameter.
#' 
#' @note The argument \code{x0} for the threshold (scale parameter) of the
#' Pareto distribution was introduced in version 0.2.
#' 
#' @author Andreas Alfons and Josef Holzer
#' 
#' @seealso \code{\link{paretoTail}}, \code{\link{fitPareto}}
#' 
#' @references Kratz, M.F. and Resnick, S.I. (1996) The QQ-estimator and heavy
#' tails.  \emph{Stochastic Models}, \bold{12}(4), 699--724.
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
#' thetaQQ(eusilc$eqIncome, k = ts$k)
#' 
#' # using threshold
#' thetaQQ(eusilc$eqIncome, x0 = ts$x0)
#' 
#' @export

thetaQQ <- function(x, k = NULL, x0 = NULL) {
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
    if(any(i <- is.na(x))) x <- x[!i]  # remove missing values
    x <- sort(x)
    n <- length(x)
#    if(haveK) {  # 'k' is supplied, threshold is determined
#        if(k >= n) stop("'k' must be smaller than the number of observed values")
#        x0 <- x[n-k]  # threshold (scale parameter)
#    } else {  # 'k' is not supplied, it is determined using threshold
#        # values are already sorted
#        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
#        k <- length(which(x > x0))
#    }
    if(!haveK) {  # 'k' is not supplied, it is determined using threshold
        # values are already sorted
        if(x0 >= x[n]) stop("'x0' must be smaller than the maximum of 'x'")
        k <- length(which(x > x0))
    }
    ## calculations
    logx <- log(x[n:(n-k+1)])  # lograrithm of reversed tail
    h <- -log((1:k)/(k+1))
    ## QQ-estimator
    (k*sum(h^2) - sum(h)^2) / sum(h * (k*logx - sum(logx)))
}
