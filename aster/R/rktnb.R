
rktnb <- function(n, size, k, mu, xpred = 1) {

    if (length(n) > 1)
        n <- length(n)
    stopifnot(n == as.integer(n))
    stopifnot(n > 0)

    stopifnot(is.numeric(xpred))
    stopifnot(length(xpred) > 0)
    stopifnot(all(xpred == as.integer(xpred)))

    stopifnot(is.numeric(mu))
    stopifnot(length(mu) > 0)
    stopifnot(all(mu > 0))

    stopifnot(is.numeric(k))
    stopifnot(length(k) >= 0)
    stopifnot(all(k == as.integer(k)))

    stopifnot(is.numeric(size))
    stopifnot(length(size) >= 0)
    stopifnot(all(size > 0))

    .C("aster_rktnb",
        n = as.integer(n),
        lenxp = length(xpred),
        lenmu = length(mu),
        lenk = length(k),
        lenalpha = length(size),
        xpred = as.double(xpred),
        mu = as.double(mu),
        k = as.integer(k),
        alpha = as.double(size),
        result = double(n), PACKAGE = "aster")$result
}

