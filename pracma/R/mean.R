##
##  m e a n . R  Geometric and Harmonic Mean (Matlab Style)
##


harmmean <- function(x, dim = 1) {
    stopifnot(is.numeric(x))
    if (dim < 1 || dim > ndims(x))
        stop("Argument 'dim' must be between 1 and 'ndims(x)'.")

    if (is.vector(x)) {
        m <- 1 / mean( 1/x )
    } else if (is.matrix(x)) {
        m <- 1 / apply( 1/x, c(2,1)[dim], mean)
    } else if (is.array(x)) {
        mid <- setdiff(1:ndims(x), dim)
        m <- 1 / apply( 1/x, mid, mean)
    } else {
        stop("Argument 'x' must be a numeric vector, matrix, or array.")
    }

    return(m)
}


geomean <- function(x, dim = 1) {
    stopifnot(is.numeric(x))
    if (dim < 1 || dim > ndims(x))
        stop("Argument 'dim' must be between 1 and 'ndims(x)'.")

    if (is.vector(x)) {
        m <- exp( sum( log(x) ) / length(x) )
    } else if (is.matrix(x)) {
        n <- size(x)[dim]
        m <- exp( apply(log(x), c(1,2)[-dim], sum) / n )
    } else if (is.array(x)) {
        n <- size(x)[dim]
        mid <- c(1:ndims(x))[-dim]
        m <- exp( apply(log(x), c(1:n)[-dim], sum) / n )
    } else {
        stop("Argument 'x' must be a numeric vector, matrix, or array.")
    }

    return(m)
}


trimmean <- function(x, percent = 0) {
    stopifnot(is.numeric(x), is.numeric(percent))
    if (length(percent) != 1 || percent < 0 || 100 < percent)
        stop("Argument 'percent' must be a scalar between 0 and 100.")

    # k <- percent / 100 / 2
    # mean(x, trim = k, na.remove = TRUE)

    .tmean <- function(x, p) {
        n <- length(x)
        k <- round (n * percent/100 / 2 )
        if (2*k > n-1) return(NA)
        x <- sort(x)
        mean(x[(k+1):(n-k)])
    }

    if (is.vector(x)) {
        m <- .tmean(x, percent)
    } else if (is.matrix(x)) {
        m <- apply(x, 2, .tmean)
    } else if (is.array(x)) {
        stop("Function 'trimmean' not yet implemented for arrays.")
    } else
        stop("Argument 'x' must be a numeric vector or matrix.")

    return(m)
}
