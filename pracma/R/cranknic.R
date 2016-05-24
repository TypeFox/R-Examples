##
##  c r a n k n i c . R  Crank-Nicolson
##


cranknic <- function(f, t0, t1, y0, ..., N = 100) {
    stopifnot(is.numeric(y0), is.numeric(t0), length(t0) == 1, 
              is.numeric(t1), length(t1) == 1)
    if (is.vector(y0)) {
        y0 <- as.matrix(y0)
    } else if (is.matrix(y0)) {
        if (ncol(y0) != 1) 
            stop("Argument 'y0' must be a row or column vector.")
    }

    fun <- match.fun(f)
    f <- function(t, y) fun(t, y, ...)

    n <- length(y0)
    y <- y0
    yout <- matrix(NA, N, n)
    yout[1, ] <- c(y0)

    h <- (t1 - t0)/(N-1)
    t <- 0
    ts <- linspace(t0, t1, N)

    # internal function used for root finding
    cnfun <- function(w)  w - y - 0.5* h * (f(t, w) + f(t, y))
    m <- length(f(t0, y0))
    if (m != n)
        stop("Function f must return a vector the same length as 'y0'.")

    # solver used for root finding
    if (n == 1) solver <- fzero
    else        solver <- fsolve

    for (i in 2:N) {
        t <- ts[i]
        w <- solver(cnfun, y)$x
        yout[i, ] <- w
        y <- w
    }

    if (n == 1) yout <- drop(yout)
    return(list(t = ts, y = yout))
}
