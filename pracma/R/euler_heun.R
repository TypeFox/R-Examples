##
##  e u l e r H e u n . R  Euler-Heun ODE Solver
##


euler_heun <- function(f, a, b, y0, n = 100, improved = TRUE, ...) {
    stopifnot(is.numeric(a), is.numeric(b), length(a) == 1, length(b) == 1,
              is.numeric(y0), length(y0) == 1)

    fun <- match.fun(f)
    f <- function(t, y) fun(t, y, ...)
    if (length(f(a, y0)) != 1)
        stop("Argument function 'f' must be an univariate function.")

    h <- (b - a)/n
    t <- seq(a, b, length.out = n+1)
    y <- numeric(n+1)
    y[1] <- y0

    for (i in 1:n) {
        y[i+1] <- y[i] + h*f(t[i], y[i])
        if (improved) {
            y[i+1] <- y[i] + h * (f(t[i], y[i]) + f(t[i+1], y[i+1]))/2.0
        }
    }
    return(list(t = t, y = y))
}
