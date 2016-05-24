##
##  p o l y A p p r o x . R  Polynomial Approximation
##


polyApprox <- function(f, a, b, n, ...) {
    if (!is.numeric(a) || !is.numeric(b) || !is.numeric(n) ||
        length(a) != 1 || length(b) != 1 || length(n) != 1 ||
        a >= b || n <= 0)
        stop("One of arguments 'a', 'b', or 'n' incorrectly chosen.")

    f1 <- match.fun(f)
    f  <- function(x) f1(x, ...)

    # Compute the Chebyshev coefficients
    cP <- chebPoly(n)
    cC <- chebCoeff(sin, a, b, n)
    p  <- drop(cC %*% cP)
    c0 <- cC[1]

    # Compute the corresponding polynomial
    q  <- c(2, -(b+a))/(b-a)
    r  <- polytrans(p, q)
    r  <- polyadd(r, c(-c0/2))

    rf <- function(x) polyval(r, x)
    ep <- fnorm(f, rf, a, b, p = Inf)
    return(list(p = p, f = rf, estim.prec = ep))
}
