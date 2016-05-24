umpu.binom <- function(x, n, p, alpha, maxiter = 10, tol = 1e-9) {

    if (! is.numeric(x)) stop("x not numeric")
    if (! is.numeric(n)) stop("n not numeric")
    if (! is.numeric(p)) stop("p not numeric")
    if (! is.numeric(alpha)) stop("alpha not numeric")
    if (! is.numeric(maxiter)) stop("maxiter not numeric")
    if (! is.numeric(tol)) stop("tol not numeric")

    if (length(n) != 1) stop("n not scalar")
    if (length(maxiter) != 1) stop("maxiter not scalar")
    if (length(tol) != 1) stop("tol not scalar")
    foo <- (length(x) > 1) + (length(p) > 1) + (length(alpha) > 1)
    if (foo > 1) stop("at most one of x, p, alpha can be non-scalar")

    if (as.integer(n) != n) stop("n not integer")
    if (as.integer(maxiter) != maxiter) stop("maxiter not integer")
    if (any(as.integer(x) != x)) stop("x not integer")

    if (! (n > 0)) stop("n not positive")
    if (! (maxiter > 0)) stop("maxiter not positive")
    if (! (tol > 0)) stop("tol not positive")
    if (! all(0 <= x & x <= n)) stop("x not in 0, ..., n")

    if (! all(0 <= p & p <= 1)) stop("p not in [0, 1]")
    if (! all(0 <= alpha & alpha <= 1)) stop("alpha not in [0, 1]")

    if (length(p) > 1) {
        out <- .C("umpubinomt",
            x = as.integer(x),
            n = as.integer(n),
            alpha = as.double(alpha),
            p = as.double(p),
            np = length(p),
            maxiter = as.integer(maxiter),
            result = double(length(p)),
            tol = as.double(tol),
            PACKAGE = "ump")
    } else if (length(alpha) > 1) {
        out <- .C("umpubinoma",
            x = as.integer(x),
            n = as.integer(n),
            alpha = as.double(alpha),
            nalpha = length(alpha),
            p = as.double(p),
            maxiter = as.integer(maxiter),
            result = double(length(alpha)),
            tol = as.double(tol),
            PACKAGE = "ump")
    } else {
        out <- .C("umpubinomx",
            x = as.integer(x),
            nx = length(x),
            n = as.integer(n),
            alpha = as.double(alpha),
            p = as.double(p),
            maxiter = as.integer(maxiter),
            result = double(length(x)),
            tol = as.double(tol),
            PACKAGE = "ump")
    }
    return(out$result)
}
