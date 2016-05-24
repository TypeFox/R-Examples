##
##  p p . R  Piecewise Polynomial Structures
##


mkpp <- function(x, P) {
    stopifnot(is.numeric(x), is.numeric(P))
    if (!is.vector(x) || is.unsorted(x))
        stop("Argument 'x' must be a sorted, finite numeric vector.")

    lx <- length(x)
    n <- nrow(P); m <- ncol(P)
    if (lx != n+1) stop("Length of 'x' must be equal to 'nrow(P)+1'.")
    
    pp <- list(breaks = x, coefs = P, pieces = n, order = m, dim = 1)
    class(pp) <- "pp"
    return(pp)
}


ppval <- function(pp, xx) {
    stopifnot(is.numeric(xx), any(!is.na(xx)))
    if (!class(pp) == "pp")
        stop("Argument 'pp' must be piecewise polynomial structure.")
    
    lx <- length(xx); yy <- rep(NA, lx)
    xb <- pp$breaks;  lb <- length(xb)

    inds <- findInterval(xx, xb)
    inds[inds == 0] <- 1
    inds[inds == lb] <- lb - 1
    for (i in 1:(lb-1)) {
        js <- which(inds == i)
        yy[js] <- polyval(pp$coefs[i, ], xx[js] - xb[i])
    }
    return(yy)
}
