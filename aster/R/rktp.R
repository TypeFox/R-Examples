
rktp <- function(n, k, mu, xpred = 1) {

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

    .C("aster_rktp",
        n = as.integer(n),
        lenxp = length(xpred),
        lenmu = length(mu),
        lenk = length(k),
        xpred = as.double(xpred),
        mu = as.double(mu),
        k = as.integer(k),
        result = double(n), PACKAGE = "aster")$result
}

rnzp <- function(n, mu, xpred = 1) rktp(n, 0, mu, xpred)

