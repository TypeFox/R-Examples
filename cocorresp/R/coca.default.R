"coca.default" <- function(y, x, method = c("predictive", "symmetric"),
                           reg.method = c("simpls", "eigen"),
                           weights = NULL,
                           n.axes = NULL,
                           symmetric = FALSE, ...) {
    nam.dat <- list(namY = deparse(substitute(y)),
                    namX = deparse(substitute(x)))
    if(any(rowSums(y) <= 0 ))
        stop("all row sums must be >0 in data matrix y")
    if(any((csum <- colSums(y)) <= 0 )) {
        y <- y[, csum > 0, drop = FALSE]
        message("some species contain no data and were removed from data matrix y\n")
    }
    if(any(rowSums(x) <= 0 ))
        stop("all row sums must be >0 in data matrix x")
    if(any((csum <- colSums(x)) <= 0 )) {
        x <- x[, csum > 0, drop = FALSE]
        message("some species contain no data and were removed from data matrix x\n")
    }
    method <- match.arg(method)
    if(method == "predictive") {
        reg.method <- match.arg(reg.method)
        retval <- switch(reg.method,
                         simpls = predcoca.simpls(y, x, R0 = weights,
                         n.axes = n.axes, nam.dat),
                         eigen = predcoca.eigen(y, x, R0 = weights,
                         n.axes = n.axes, nam.dat))
    } else {
        retval <- symcoca(y, x, n.axes = n.axes, R0 = weights,
                          symmetric = symmetric, nam.dat)
    }
    class(retval) <- c("coca", class(retval))
    retval
}

