##
##  p p f i t . R  Piecewise Polynomial Fit
##


ppfit <- function(x, y, xi, method = c("linear", "cubic")) {
    stopifnot(is.numeric(x), is.numeric(y), is.numeric(xi))
    if (length(x) != length(y) || length(x) <= 1)
        stop("Length of 'x' and 'y' must be equal and greater than 1.")
    method <- match.arg(method)

    y0 <- interp1(x, y, xi, method = method)
    nn <- finds(x >= xi[1] & x <= xi[length(xi)])
    fcn <- function(yi) sum((interp1 (xi, yi, x[nn], method) - y[nn])^2)
    yi <- optim(y0, fcn)$par

    n <- length(xi) - 1
    if (method == "linear") {
        P <- matrix(NA, nrow = n, ncol = 2)
        for (i in 1:n)
            P[i, ] <- c((yi[i+1]-yi[i])/(xi[i+1]-xi[i]), yi[i])
        pp <- mkpp(xi, P)
    } else if (method == "cubic") {
        pp <- cubicspline(xi, yi)
    } else
        stop("Unknown method: must be 'linear' or 'cubic'.")

    return(pp)
}
