##
##  f n o r m .R  Function Norm
##


fnorm <- function(f, g, x1, x2, p = 2, npoints = 100) {
    stopifnot(is.numeric(x1), length(x1) == 1,
              is.numeric(x2), length(x2) == 1, x1 < x2,
              is.numeric(npoints), length(npoints) == 1, npoints >= 2)

    f <- match.fun(f)
    g <- match.fun(g)

    x  <- seq(x1, x2, length.out = npoints)
    yf <- f(x)
    yg <- g(x)
    if (length(yf) != npoints || length(yg) != npoints)
        stop("Arguments 'f' and 'g' must be vectorized functions.")

    fd <- Norm(yf - yg, p = p)
    return(fd)
}
