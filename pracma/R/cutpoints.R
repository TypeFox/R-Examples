cutpoints <- function(x, nmax = 8, quant = 0.95) {
    stopifnot(is.numeric(x), is.numeric(nmax), is.numeric(quant))
    if (length(nmax) != 1 || floor(nmax) != ceiling(nmax) || nmax <= 0)
        stop("Argument 'nmax' must be a positive integer.")
    if (length(quant) != 1 || quant < 0 || quant > 1)
        stop("Argument 'quant' must be a numeric value in [0,1].")

    if (length(x) == 1) return(list(cutp = c(), cutd = c()))
    if (!is.vector(x))  x <- c(x)
    if (is.unsorted(x)) x <- sort(x)

    d <- diff(x)
    if (quant == 0.0) n <- nmax
    else {
        q <- quantile(d, probs = quant)
        n <- min(nmax, sum(d >= q))
    }

    o <- order(d)
    inds <- o[length(d):(length(d)-(n-1))]
    dc <- d[inds]
    xc <- x[inds] + dc/2

    o <- order(xc)
    return(list(cutp = xc[o], cutd = dc[o]))
}

##  Example
