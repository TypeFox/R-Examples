
# ---------------------------------------
# Author: Stefan Kraft
#         Vienna University of Technology
# ---------------------------------------

# Code copied from "simPopulation" package version 0.4.0

quantileWt <- function(x, weights = NULL, 
        probs = seq(0, 1, 0.25), na.rm = TRUE) {
    # initializations
    if(!is.numeric(x)) stop("'x' must be a numeric vector")
    x <- unname(x)  # unlike 'quantile', this never returns a named vector
    if(is.null(weights)) {
        return(quantile(x, probs, na.rm=na.rm, names=FALSE, type=1))
    } else if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
    else if(length(weights) != length(x)) {
        stop("'weights' must have the same length as 'x'")
    } else if(!all(is.finite(weights))) stop("missing or infinite weights")
    if(!is.numeric(probs) || all(is.na(probs)) || 
            isTRUE(any(probs < 0 | probs > 1))) {
        stop("'probs' must be a numeric vector with values in [0,1]")
    }
    if(length(x) == 0) return(rep.int(NA, length(probs)))
    if(!isTRUE(na.rm) && any(is.na(x))) {
        stop("missing values and NaN's not allowed if 'na.rm' is not TRUE")
    }
    # sort values and weights
    ord <- order(x, na.last=NA)
    x <- x[ord]
    weights <- weights[ord]
    # some preparations
    rw <- cumsum(weights)/sum(weights)
    # obtain quantiles
    select <- sapply(probs, function(p) min(which(rw >= p)))
    q <- x[select]
    return(q)
}

spBwplotStats <- function(x, weights = NULL, coef = 1.5, 
        zeros = TRUE, do.out = TRUE) {
    # initializations
    if(!is.numeric(x)) stop("'x' must be a numeric vector")
    if(!is.numeric(coef) || length(coef) != 1 || coef < 0) {
        stop("'coef' must be a single non-negative number")
    }
    # get quantiles
    if(isTRUE(zeros)) {
        zero <- ifelse(is.na(x), FALSE, x == 0)
        x <- x[!zero]
        if(is.null(weights)) nzero <- sum(zero)
        else {
            # if 'zeros' is not TRUE, these checks are done in 'quantileWt'
            # but here we need them since we use subscripting
            if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
            else if(length(weights) != length(zero)) {
                stop("'weights' must have the same length as 'x'")
            }
            nzero <- sum(weights[zero])
            weights <- weights[!zero]
        }
    } else nzero <- NULL
    ok <- !is.na(x)
    n <- if(is.null(weights)) sum(ok) else sum(weights[ok])
    if(n == 0) stats <- rep.int(NA, 5)
    else stats <- quantileWt(x, weights)
    iqr <- diff(stats[c(2, 4)])  # inter quartile range
    if(coef == 0) do.out <- FALSE
    else {
        if(is.na(iqr)) out <- is.infinite(x) 
        else {
            lower <- stats[2] - coef * iqr
            upper <- stats[4] + coef * iqr
            out <- ifelse(ok, x < lower | x > upper, FALSE)
        }
        if(any(out)) stats[c(1, 5)] <- range(x[!out], na.rm=TRUE)
    }
    res <- list(stats=stats, n=n, nzero=nzero, 
        out=if(isTRUE(do.out)) x[out] else numeric())
    class(res) <- "spBwplotStats"
    res
}
