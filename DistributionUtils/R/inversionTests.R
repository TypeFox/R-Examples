inversionTestpq <- function(densFn = "norm", n = 10,
                            intTol = .Machine$double.eps^0.25,
                            uniTol = intTol, x = NULL, method = "spline", ...)
{
    if (is.null(x)){
        rfun <- match.fun(paste("r", densFn, sep = ""))
        x <- rfun(n, ...)
    } else x <- x
    qpx <- qDist(densFn, p = pDist(densFn, q = x, intTol = intTol, ...),
                 uniTol = uniTol, method = "spline", ...)
    diffs <- qpx - x
    result <- list(qpx = qpx, x = x, diffs = diffs, n = length(x))

    return(result)
}

inversionTestqp <- function(densFn = "norm",
                            qs = c(0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4,
                                   0.5, 0.6, 0.8, 0.9, 0.95, 0.975, 0.99,
                                   0.999),
                            uniTol = .Machine$double.eps^0.25,
                            intTol = uniTol, method = "spline", ...)
{
    q <- qDist(densFn, p = qs, uniTol = uniTol, ...)
    pqqs <- pDist(densFn, q = q, intTol = intTol, ...)
    diffs <- pqqs - qs
    result <- list(pqqs = pqqs, qs = qs, diffs = diffs)

    return(result)
}
