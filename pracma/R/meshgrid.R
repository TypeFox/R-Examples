##
##  m e s h g r i d . R  Generate a Mesh Grid
##


meshgrid <- function(x, y = x) {
    if (!is.numeric(x) || !is.numeric(y))
        stop("Arguments 'x' and 'y' must be numeric vectors.")

    x <- c(x); y <- c(y)
    n <- length(x)
    m <- length(y)

    X <- matrix(rep(x, each = m),  nrow = m, ncol = n)
    Y <- matrix(rep(y, times = n), nrow = m, ncol = n)

    return(list(X = X, Y = Y))
}


peaks <- function(v = 49, w) {
    stopifnot(is.numeric(v))
    if (missing(w)) {
        if (length(v) == 1 && v >= 1) {
            mg <- meshgrid(linspace(-3, 3, floor(v)))
             x <- mg$X; y <- mg$Y
        } else {
            mg <- meshgrid(v, v)
             x <- mg$X; y <- mg$Y
        }
    } else {
        stopifnot(is.numeric(w))
        x <- v; y <- w
    }

    z <-  3 * (1-x)^2 * exp(-(x^2) - (y+1)^2) -
         10 * (x/5 - x^3 - y^5) * exp(-x^2 - y^2) -
        1/3 * exp(-(x+1)^2 - y^2)

    return(list(X = x, Y = y, Z = z))
}
