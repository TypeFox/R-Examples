# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

#' Fit income distribution models with the Pareto distribution
#' 
#' Fit a Pareto distribution to the upper tail of income data.  Since a
#' theoretical distribution is used for the upper tail, this is a semiparametric
#' approach.
#' 
#' The arguments \code{k} and \code{x0} of course correspond with each other.
#' If \code{k} is supplied, the threshold \code{x0} is estimated with the \eqn{n
#' - k} largest value in \code{x}, where \eqn{n} is the number of observations.
#' On the other hand, if the threshold \code{x0} is supplied, \code{k} is given
#' by the number of observations in \code{x} larger than \code{x0}.  Therefore,
#' either \code{k} or \code{x0} needs to be supplied.  If both are supplied,
#' only \code{k} is used (mainly for back compatibility).
#' 
#' The function supplied to \code{method} should take a numeric vector (the
#' observations) as its first argument.  If \code{k} is supplied, it will be
#' passed on (in this case, the function is required to have an argument called
#' \code{k}).  Similarly, if the threshold \code{x0} is supplied, it will be
#' passed on (in this case, the function is required to have an argument called
#' \code{x0}).  As above, only \code{k} is passed on if both are supplied.  If
#' the function specified by \code{method} can handle sample weights, the
#' corresponding argument should be called \code{w}.  Additional arguments are
#' passed via the \dots{} argument.
#' 
#' @param x a numeric vector.
#' @param k the number of observations in the upper tail to which the Pareto
#' distribution is fitted.
#' @param x0 the threshold (scale parameter) above which the Pareto distribution
#' is fitted.
#' @param method either a function or a character string specifying the function
#' to be used to estimate the shape parameter of the Pareto distibution, such as
#' \code{\link{thetaPDC}} (the default).  See \dQuote{Details} for requirements
#' for such a function and \dQuote{See also} for available functions.
#' @param groups an optional vector or factor specifying groups of elements of
#' \code{x} (e.g., households).  If supplied, each group of observations is
#' expected to have the same value in \code{x} (e.g., household income).  Only
#' the values of every first group member to appear are used for fitting the
#' Pareto distribution. For each group above the threshold, every group member
#' is assigned the same value.
#' @param w an optional numeric vector giving sample weights.
#' @param \dots addtional arguments to be passed to the specified method.
#' 
#' @return A numeric vector with a Pareto distribution fit to the upper tail.
#' 
#' @note The arguments \code{x0} for the threshold (scale parameter) of the
#' Pareto distribution and \code{w} for sample weights were introduced in
#' version 0.2.  This results in slightly different behavior regarding the
#' function calls to \code{method} compared to prior versions.
#' 
#' @author Andreas Alfons and Josef Holzer
#' 
#' @seealso \code{\link{paretoTail}}, \code{\link{replaceTail}}
#' 
#' \code{\link{thetaPDC}}, \code{\link{thetaWML}}, \code{\link{thetaHill}},
#' \code{\link{thetaISE}}, \code{\link{thetaLS}}, \code{\link{thetaMoment}},
#' \code{\link{thetaQQ}}, \code{\link{thetaTM}}
#' 
#' @keywords manip
#' 
#' @examples
#' data(eusilc)
#' 
#' 
#' ## gini coefficient without Pareto tail modeling
#' gini("eqIncome", weights = "rb050", data = eusilc)
#' 
#' 
#' ## gini coefficient with Pareto tail modeling
#' 
#' # using number of observations in tail
#' eqIncome <- fitPareto(eusilc$eqIncome, k = 175, 
#'     w = eusilc$db090, groups = eusilc$db030)
#' gini(eqIncome, weights = eusilc$rb050)
#' 
#' # using threshold
#' eqIncome <- fitPareto(eusilc$eqIncome, x0 = 44150, 
#'     w = eusilc$db090, groups = eusilc$db030)
#' gini(eqIncome, weights = eusilc$rb050)
#' 
#' @export

fitPareto <- function(x, k = NULL, x0 = NULL, 
        method = "thetaPDC", groups = NULL, w = NULL, ...) {
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
    nam <- argNames(method)
    useW <- !is.null(w) && ("w" %in% nam)
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
        values <- x[unique]
        if(useW) w <- w[unique]
    } else values <- x
    ## check for missing values
    indices <- 1:length(values)
    if(any(i <- is.na(values))) indices <- indices[!i]
    ## order of observed values
    order <- order(values[indices])
    indicesSorted <- indices[order]  # indices of sorted vector
    n <- length(indicesSorted)
    ## start constructing call to 'method' for estimation of shape parameter
    dots <- list(values[indices], ...)
    if(haveK) {  # 'k' is supplied, threshold is determined
        if(k >= n) stop("'k' must be smaller than the number of observed values")
        x0 <- values[indicesSorted[n-k]]  # threshold (scale parameter)
        dots$k <- k  # 'method' is expected to have 'k' as argument
    } else {  # 'k' is not supplied, it is determined using threshold
        if(x0 >= values[indicesSorted[n]]) {  # compare to sorted values
            stop("'x0' must be smaller than the largest value")
        }
        k <- length(which(values[indices] > x0))  # number of observations in tail
        dots$x0 <- x0  # 'method' is expected to have threshold 'x0' as argument
    }
    ## estimate shape parameter
    if(useW) dots$w <- w[indices]
    theta <- do.call(method, dots)
    ## fit Pareto distribution
    valuesPareto <- x0/runif(k)^(1/theta)
    values[indicesSorted[(n-k+1):n]] <- sort(valuesPareto)
    ## return values
    if(haveGroups) {
        groups <- as.character(groups)
        names(values) <- groups[unique]
        values <- values[groups]
        names(values) <- names(x)
    }
    values
}
