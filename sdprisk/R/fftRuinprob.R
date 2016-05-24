fftRuinprob <- function(process, interval, maxreserve, n, use.splines = FALSE) {
    stopifnot(is.logical(use.splines))

    misspar <- c(missing(interval),
                 missing(maxreserve),
                 missing(n))

    if (sum(as.integer(misspar)) > 1L) {
        stop('At least two parameters of ', sQuote('interval'), ', ',
             sQuote('maxreserve'), ' or ', sQuote('n'), ' are required.')
    } else {
        if (misspar[3L]) {
            n <- nextn(maxreserve / interval)
        }

        if (misspar[2L]) {
            n          <- nextn(n)
            maxreserve <- n * interval
        }

        if (misspar[1L]) {
            n        <- nextn(n)
            interval <- maxreserve / n
        }

        stopifnot(is.numeric(interval),
                  is.numeric(maxreserve),
                  is.numeric(n))
    }

    p      <- process[['p']]
    q      <- process[['q']]
    zeta   <- process[['zeta']]

    x <- seq.int(from       = 0.0,
                  by         = interval,
                  length.out = n + 1L)

    if (is.finite(zeta)) {
        a.fft <- fft(diff(pexp(q    = x,
                               rate = zeta)))
    } else {
        a.fft <- rep.int(1.0 + 0.0i, n)
    }

    b.fft <- fft(diff(vapply(X         = x,
                             FUN       = process[[c('claims', 'cdf.tailarea')]],
                             FUN.VALUE = numeric(1L))))

    w <- Re(fft(z       = p * a.fft / (1 - q * a.fft * b.fft),
                inverse = TRUE) / n)

    x  <- x[-(n + 1L)]
    y  <- rev(cumsum(rev(w)))

    if (!is.finite(zeta)) {
        y[1L] <- q
        y1    <- rep.int(0.0, n)
    } else {
        y1    <- w / (zeta * p * interval)
    }

    if (use.splines) {
        psi   <- splinefun(x, y)
        psi.1 <- splinefun(x, y1)
    } else {
        psi   <- approxfun(x, y,  rule = 2L)
        psi.1 <- approxfun(x, y1, rule = 2L)
    }

    return(structure(.Data       = list(psi   = psi,
                                        psi.1 = psi.1,
                                        psi.2 = function(x) psi(x) - psi.1(x)),
                     compmethod  = 'fft',
                     riskproc    = process,
                     parameters  = list(interval    = interval,
                                        maxreserve  = maxreserve,
                                        n           = n,
                                        use.splines = use.splines),
                     diagnostics = list(w = w)))
}
