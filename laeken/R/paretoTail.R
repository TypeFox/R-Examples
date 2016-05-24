# ----------------------------------------
# Authors: Andreas Alfons and Josef Holzer
#          Vienna University of Technology
# ----------------------------------------

#' Pareto tail modeling for income distributions
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
#' only \code{k} is used.
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
#' @aliases print.paretoTail
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
#' Pareto distribution.
#' @param w an optional numeric vector giving sample weights.
#' @param alpha numeric; values above the theoretical \eqn{1 - }\code{alpha}
#' quantile of the fitted Pareto distribution will be flagged as outliers for
#' further treatment with \code{\link{reweightOut}} or \code{\link{replaceOut}}.
#' @param \dots addtional arguments to be passed to the specified method.
#'
#' @returnClass paretoTail
#' @returnItem x the supplied numeric vector.
#' @returnItem k the number of observations in the upper tail to which the
#' Pareto distribution has been fitted.
#' @returnItem groups if supplied, the vector or factor specifying groups of
#' elements.
#' @returnItem w if supplied, the numeric vector of sample weights.
#' @returnItem method the function used to estimate the shape parameter, or the
#' name of the function.
#' @returnItem x0 the scale parameter.
#' @returnItem theta the estimated shape parameter.
#' @returnItem tail if \code{groups} is not \code{NULL}, this gives the groups
#' with values larger than the threshold (scale parameter), otherwise the
#' indices of observations in the upper tail.
#' @returnItem alpha the tuning parameter \code{alpha} used for flagging
#' outliers.
#' @returnItem out if \code{groups} is not \code{NULL}, this gives the groups
#' that are flagged as outliers, otherwise the indices of the flagged
#' observations.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{reweightOut}}, \code{\link{shrinkOut}},
#' \code{\link{replaceOut}}, \code{\link{replaceTail}}, \code{\link{fitPareto}}
#'
#' \code{\link{thetaPDC}}, \code{\link{thetaWML}}, \code{\link{thetaHill}},
#' \code{\link{thetaISE}}, \code{\link{thetaLS}}, \code{\link{thetaMoment}},
#' \code{\link{thetaQQ}}, \code{\link{thetaTM}}
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
#' # estimate threshold
#' ts <- paretoScale(eusilc$eqIncome, w = eusilc$db090,
#'     groups = eusilc$db030)
#'
#' # estimate shape parameter
#' fit <- paretoTail(eusilc$eqIncome, k = ts$k,
#'     w = eusilc$db090, groups = eusilc$db030)
#'
#' # calibration of outliers
#' w <- reweightOut(fit, calibVars(eusilc$db040))
#' gini(eusilc$eqIncome, w)
#'
#' # winsorization of outliers
#' eqIncome <- shrinkOut(fit)
#' gini(eqIncome, weights = eusilc$rb050)
#'
#' # replacement of outliers
#' eqIncome <- replaceOut(fit)
#' gini(eqIncome, weights = eusilc$rb050)
#'
#' # replacement of whole tail
#' eqIncome <- replaceTail(fit)
#' gini(eqIncome, weights = eusilc$rb050)
#'
#' @export

paretoTail <- function(x, k = NULL, x0 = NULL, method = "thetaPDC",
        groups = NULL, w = NULL, alpha = 0.01, ...) {
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
    if(is.character(method)) method <- getDotTheta(method)
    nam <- argNames(method)
    useW <- !is.null(w) && ("w" %in% nam)
    if(useW && (!is.numeric(w) || length(w) != length(x))) {
        stop("'w' must be numeric vector of the same length as 'x'")
    }
    ## allow for cluster effect
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
        xx <- x[unique]
        if(useW) ww <- w[unique]
    } else {
        xx <- x
        if(useW) ww <- w
    }
    xx <- unname(xx)
    ## check for missing values
    if(any(i <- is.na(xx))) {
        xx <- xx[!i]
        if(useW) ww <- ww[!i]
    }
    ## order of observed values
    order <- order(xx)
    xx <- xx[order]
    if(useW) ww <- ww[order]
    n <- length(xx)
    ## start constructing call to 'method' for estimation of shape parameter
    dots <- list(xx, ...)
    if(haveK) {  # 'k' is supplied, threshold is determined
        if(k >= n) stop("'k' must be smaller than the number of observed values")
        x0 <- xx[n-k]  # threshold (scale parameter)
        dots$k <- k  # 'method' is expected to have 'k' as argument
    } else {  # 'k' is not supplied, it is determined using threshold
        if(x0 >= xx[n]) {  # compare to sorted values
            stop("'x0' must be smaller than the largest value")
        }
        k <- length(which(xx > x0))  # number of observations in tail
        dots$x0 <- x0  # 'method' is expected to have threshold 'x0' as argument
    }
    ## estimate shape parameter
    if(useW) dots$w <- ww
    theta <- do.call(method, dots)
    ## indicate observations in tail
    if(haveGroups) {
        tail <- groups[unique]
        tail <- tail[!i]
        tail <- tail[order]
        tail <- tail[(n-k+1):n]
    } else tail <- order[(n-k+1):n]
    ## flag suspicious observations (nonrepresentative outliers)
    if(!is.numeric(alpha) || length(alpha) == 0 || alpha < 0 || alpha > 1) {
        stop("'alpha' must be a numeric value in [0,1]")
    } else alpha <- alpha[1]
    q <- qpareto(1-alpha, x0, theta)  # quantile of the Pareto distribution
    if(haveGroups) {
        out <- which(xx[(n-k+1):n] > q)
        out <- tail[out]
    } else {
        out <- unname(which(x > q))
        out <- out[order(x[out])]
    }

    ## return object
    res <- list(x=x, k=k, groups=groups, w=w, method=method,
        x0=x0, theta=theta, tail=tail, alpha=alpha, out=out)
    class(res) <- "paretoTail"
    res
}


#' Replace observations under a Pareto model
#'
#' Replace observations under a Pareto model for the upper tail with values
#' drawn from the fitted distribution.
#'
#' \code{replaceOut(x, \dots{})} is a simple wrapper for \code{replaceTail(x,
#' all = FALSE, \dots{})}.
#'
#' @param x an object of class \code{"paretoTail"} (see
#' \code{\link{paretoTail}}).
#' @param all a logical indicating whether all observations in the upper tail
#' should be replaced or only those flagged as outliers.
#' @param \dots additional arguments to be passed down.
#'
#' @return A numeric vector consisting mostly of the original values, but with
#' observations in the upper tail replaced with values from the fitted Pareto
#' distribution.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{paretoTail}}, \code{\link{reweightOut}},
#' \code{\link{shrinkOut}}
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
#' # estimate threshold
#' ts <- paretoScale(eusilc$eqIncome, w = eusilc$db090,
#'     groups = eusilc$db030)
#'
#' # estimate shape parameter
#' fit <- paretoTail(eusilc$eqIncome, k = ts$k,
#'     w = eusilc$db090, groups = eusilc$db030)
#'
#' # replacement of outliers
#' eqIncome <- replaceOut(fit)
#' gini(eqIncome, weights = eusilc$rb050)
#'
#' # replacement of whole tail
#' eqIncome <- replaceTail(fit)
#' gini(eqIncome, weights = eusilc$rb050)
#'
#' @export

replaceTail <- function(x, ...) UseMethod("replaceTail")

#' @rdname replaceTail
#' @method replaceTail paretoTail
#' @export
replaceTail.paretoTail <- function(x, all = TRUE, ...) {
    which <- if(isTRUE(all)) x$tail else x$out
    k <- length(which)  # number of observations to be replaced
    res <- x$x
    if(k > 0) {
        new <- sort(rpareto(k, x$x0, x$theta))
        groups <- x$groups
        if(is.null(groups)) res[which] <- new
        else {
            groups <- as.character(groups)
            which <- as.character(which)
            replace <- which(groups %in% which)
            names(new) <- which
            new <- new[groups[replace]]
            names(new) <- names(res[replace])
            res[replace] <- new
        }
    }
    res
}

#replaceOut <- function(x) UseMethod("replaceOut")
#
#replaceOut.paretoTail <- function(x) {
#    out <- x$out
#    nout <- length(out)  # number of nonrepresentative outliers
#    new <- sort(rpareto(nout, x$x0, x$theta))
#    res <- x$x
#    if(nout > 0) {
#        groups <- x$groups
#        if(is.null(groups)) res[out] <- new
#        else {
#            groups <- as.character(groups)
#            out <- as.character(out)
#            replace <- which(groups %in% out)
#            names(new) <- out
#            new <- new[groups[replace]]
#            names(new) <- names(res[replace])
#            res[replace] <- new
#        }
#    }
#    res
#}

#replaceOut <- function(x) replaceTail(x, all=FALSE)

#' @rdname replaceTail
#' @export
replaceOut <- function(x, ...) {
    localReplaceTail <- function(x, all, ...) replaceTail(x, all=FALSE, ...)
    localReplaceTail(x, ...)
}


#' Reweight outliers in the Pareto model
#'
#' Reweight observations that are flagged as outliers in a Pareto model for the
#' upper tail of the distribution.
#'
#' If the data contain sample weights, the weights of the outlying observations
#' are set to \eqn{1} and the weights of the remaining observations are
#' calibrated according to auxiliary variables.  Otherwise, weight \eqn{0} is
#' assigned to outliers and weight \eqn{1} to other observations.
#'
#' @param x an object of class \code{"paretoTail"} (see
#' \code{\link{paretoTail}}).
#' @param X a matrix of binary calibration variables (see
#' \code{\link{calibVars}}).  This is only used if \code{x} contains sample
#' weights or if \code{w} is supplied.
#' @param w a numeric vector of sample weights. This is only used if \code{x}
#' does not contain sample weights, i.e., if sample weights were not considered
#' in estimating the shape parameter of the Pareto distribution.
#' @param \dots additional arguments to be passed down.
#'
#' @return If the data contain sample weights, a numeric containing the
#' recalibrated weights is returned, otherwise a numeric vector assigning weight
#' \eqn{0} to outliers and weight \eqn{1} to other observations.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{paretoTail}}, \code{\link{shrinkOut}} ,
#' \code{\link{replaceOut}}, \code{\link{replaceTail}}
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
#' @keywords manip
#'
#' @examples
#' data(eusilc)
#'
#' ## gini coefficient without Pareto tail modeling
#' gini("eqIncome", weights = "rb050", data = eusilc)
#'
#' ## gini coefficient with Pareto tail modeling
#' # estimate threshold
#' ts <- paretoScale(eusilc$eqIncome, w = eusilc$db090,
#'     groups = eusilc$db030)
#' # estimate shape parameter
#' fit <- paretoTail(eusilc$eqIncome, k = ts$k,
#'     w = eusilc$db090, groups = eusilc$db030)
#' # calibration of outliers
#' w <- reweightOut(fit, calibVars(eusilc$db040))
#' gini(eusilc$eqIncome, w)
#'
#' @export

reweightOut <- function(x, ...) UseMethod("reweightOut")

#' @rdname reweightOut
#' @method reweightOut paretoTail
#' @export
reweightOut.paretoTail <- function(x, X, w = NULL, ...) {
    # in case of sample weights, set weights of outliers to one and calibrate
    # other observations
    # otherwise, set weights of outliers to zero and weights of other
    # observations to one
    out <- x$out
    n <- length(x$x)  # number of observations
    if(is.null(x$w)) {
        # check supplied weights
        if(!is.null(w) && (!is.numeric(w) || length(w) != n)) {
            stop(sprintf("'w' must be numeric vector of length %d", n))
        }
    } else w <- x$w
    if(length(out) > 0) {  # nonrepresentative outliers
        groups <- x$groups
        if(!is.null(groups)) out <- which(groups %in% out)
        if(is.null(w)) {
            w <- rep.int(1, n)
            w[out] <- 0
        } else {
            totals <- apply(X, 2, function(i) sum(i*w))
            args <- list(...)
            args$X <- X[-out, , drop=FALSE]
            args$d <- w[-out]
            w[out] <- 1  # set weight of nonrepresentative outliers to 1
            totalsOut <- apply(X[out, , drop=FALSE], 2, sum)
            args$totals <- totals - totalsOut
            g <- do.call("calibWeights", args)
            w[-out] <- g * args$d
        }
    }
    w
}


#' Shrink outliers in the Pareto model
#'
#' Shrink observations that are flagged as outliers in a Pareto model for the
#' upper tail of the distribution to the theoretical quantile used for outlier
#' detection.
#'
#' @param x an object of class \code{"paretoTail"} (see
#' \code{\link{paretoTail}}).
#' @param \dots additional arguments to be passed down (currently ignored as
#' there are no additional arguments in the only method implemented).
#' @return A numeric vector consisting mostly of the original values, but with
#' outlying observations in the upper tail shrunken to the corresponding
#' theoretical quantile of the fitted Pareto distribution.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{paretoTail}}, \code{\link{reweightOut}},
#' \code{\link{replaceOut}}, \code{\link{replaceTail}}
#'
#' @references
#' A. Alfons and M. Templ (2013) Estimation of Social Exclusion Indicators
#' from Complex Surveys: The \R Package \pkg{laeken}.  \emph{Journal of
#' Statistical Software}, \bold{54}(15), 1--25.  URL
#' \url{http://www.jstatsoft.org/v54/i15/}
#'
#' @keywords manip
#'
#' @examples
#' data(eusilc)
#'
#' ## gini coefficient without Pareto tail modeling
#' gini("eqIncome", weights = "rb050", data = eusilc)
#'
#' ## gini coefficient with Pareto tail modeling
#' # estimate threshold
#' ts <- paretoScale(eusilc$eqIncome, w = eusilc$db090,
#'     groups = eusilc$db030)
#' # estimate shape parameter
#' fit <- paretoTail(eusilc$eqIncome, k = ts$k,
#'     w = eusilc$db090, groups = eusilc$db030)
#' # shrink outliers
#' eqIncome <- shrinkOut(fit)
#' gini(eqIncome, weights = eusilc$rb050)
#'
#' @export

shrinkOut <- function(x, ...) UseMethod("shrinkOut")

#' @rdname shrinkOut
#' @method shrinkOut paretoTail
#' @export
shrinkOut.paretoTail <- function(x, ...) {
    # winsorize outliers in the upper tail
    out <- x$out
    res <- x$x
    if(length(out) > 0) {  # nonrepresentative outliers
        new <- qpareto(1-x$alpha, x$x0, x$theta)  # quantile of the Pareto distribution
        groups <- x$groups
        if(!is.null(groups)) out <- which(groups %in% out)
        res[out] <- new
    }
    res
}


## print method for class "paretoTail"
#' @export
print.paretoTail <- function(x, ...) {
    cat("Threshold: ")
    cat(x$x0, ...)
    items <- if(is.null(x$groups)) "observations" else "groups"
    cat(sprintf("\nNumber of %s in the tail: ", items))
    cat(x$k, ...)
    cat("\nShape parameter: ")
    cat(x$theta, ...)
    cat(sprintf("\n\nOutlying %s:\n", items))
    print(x$out, ...)
}


## utility functions for Pareto distribution
dpareto <- function(x, x0 = 1, theta = 1) theta*x0^theta / x^(theta+1)
ppareto <- function(q, x0 = 1, theta = 1) 1 - (q/x0)^(-theta)
qpareto <- function(p, x0 = 1, theta = 1) unname(x0*exp(qexp(p)/theta))
rpareto <- function(n, x0 = 1, theta = 1) x0/runif(n)^(1/theta)


## other utility functions
getDotTheta <- function(method) {
    if(length(method) == 0) stop("'method' has length 0")
    else method <- method[1]
    if(method %in% c("thetaPDC", "thetaISE", "thetaWML", "thetaHill")) {
        method <- paste(".", method, sep="")
    }
    method
}
