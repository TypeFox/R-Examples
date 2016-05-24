##
##  c u r v e f i t . R  Polynomial Curve Fit
##


curvefit <- function(u, x, y, n, U = NULL, V = NULL) {
    stopifnot(is.numeric(u), is.numeric(x), is.numeric(y))
    u <- c(u); m <- length(u)
    x <- as.matrix(c(x)); y <- as.matrix(c(y))
    if (length(x) != m || length(y) != m)
        stop("Vectors 'x' and 'y' must have the same length as 't'.")
    if (!is.numeric(n) || length(n) != 1 || floor(n) != ceiling(n) || n < 1)
        stop("Argument 'n' must be an integer greater or equal 1.")

    Ax <- outer(u, seq(n, 0), "^")
    Ay <- outer(u, seq(n, 0), "^")

    if (is.null(U) || is.null(V)) {
        px <- lsqlin(Ax, x)
        py <- lsqlin(Ay, y)
    } else {
        Cx <- outer(U, seq(n, 0), "^"); cx <- V[, 1]
        Cy <- outer(U, seq(n, 0), "^"); cy <- V[, 2]
        px <- lsqlin(Ax, x, Cx, cx)
        py <- lsqlin(Ay, y, Cy, cy)
    }
    xp <- polyval(c(px), u)
    yp <- polyval(c(py), u)

    return(list(xp = xp, yp = yp, px = px, py = py))
}
