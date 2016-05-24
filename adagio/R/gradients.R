##
##  g r a d i e n t . R  Numerical Gradient and Jacobian
##


ns.grad <-
function (fn, x0, dir = c("central", "forward", "backward"),
          eps = .Machine$double.eps, ...) 
{
    if (!is.numeric(x0)) 
        stop("Argument 'x0' must be a numeric vector.")

    fct <- match.fun(fn)
    fn  <- function(x) fct(x, ...)
    if (length(fn(x0)) != 1) 
        stop("Function 'f' must be a univariate function of 2 variables.")

    dir <- match.arg(dir)
    if (dir == "central") heps <- eps^(1/3)
    else                  heps <- eps^(1/2)

    n <- length(x0)
    hh <- rep(0, n)
    gr <- numeric(n)
    for (i in 1:n) {
        hh[i] <- heps
        if (dir == "forward") {
            gr[i] <- (fn(x0 + hh) - fn(x0)) / heps
        } else if (dir == "backward") {
            gr[i] <- (fn(x0) - fn(x0 - hh)) / heps
        } else if (dir == "central") {
            gr[i] <- (fn(x0 + hh) - fn(x0 - hh)) / (2*heps)
        } else {
            stop("Direction must be one of 'forward', 'backward', or 'central'.")
        }
        hh[i] <- 0
    }

    gr
}
