##
##  q u a d c c . R  Adaptive Clenshaw-Curtis Quadrature
##


quadcc <- function(f, a, b, tol = .Machine$double.eps^0.5, ...) {
    stopifnot(is.numeric(a), length(a) == 1,
              is.numeric(b), length(b) == 1)
    eps <- .Machine$double.eps

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)
    
    if (a == b)     return(0)
    else if (a > b) return(-1 * quadgk(f, b, a, tol = tol))
    if (!is.finite(f(a))) a <- a + eps * sign(b-a)
    if (!is.finite(f(b))) b <- b - eps * sign(b-a)

    .ccadpt <- function(f, a, b, tol = tol) {
        Q4 <- clenshaw_curtis(f, a, b, n = 4)
        Q8 <- clenshaw_curtis(f, a, b, n = 8)
        if (abs(Q4 - Q8) < tol) {
            return(Q8)
        } # else

        Q2 <- .ccadpt(f, (a+b)/2, b, tol = tol)
        Q1 <- .ccadpt(f, a, (a+b)/2, tol = tol)

        return(Q1 + Q2)
    }
    return(.ccadpt(f, a, b, tol))
}

