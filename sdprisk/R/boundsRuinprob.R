boundsRuinprob <- function(process, interval, maxreserve, richardson = TRUE, use.splines = FALSE) {
    stopifnot(is.riskproc(process),
              is.numeric(interval),
              is.numeric(maxreserve),
              is.logical(richardson),
              is.logical(use.splines))

    p      <- process[['p']]
    q      <- process[['q']]
    zeta   <- process[['zeta']]
    claims <- process[['claims']]

    n <- floor(maxreserve / interval) + 3L

    mytailarea <- claims[['cdf.tailarea']]

    if (!is.null(mytailarea)) {
        h2l <- diff(vapply(X         = seq.int(from       = 0.0,
                                               by         = interval,
                                               length.out = n + 1L),
                           FUN       = mytailarea,
                           FUN.VALUE = numeric(1L)))
    } else {
        warning('Integrated tail area distribution function not supplied.\n',
                'Trying to use a rough approximation based on the CDF.\n',
                immediate. = TRUE,
                call.      = FALSE)

        mu    <- mean(claims)
        mycdf <- claims[['cdf']]

        if (is.finite(mu) & !is.null(mycdf)) {
            h2l <- (1.0 - vapply(X         = seq.int(from       = interval,
                                                     by         = interval,
                                                     length.out = n),
                                 FUN       = mycdf,
                                 FUN.VALUE = numeric(1L))) * interval / mu
        } else {
            stop('Claim CDF or mean claim size not available. Please refer to the help.',
                 call. = FALSE)
        }
    }

    if (is.finite(zeta)) {
        h1l <- diff(pexp(q    = seq.int(from       = 0.0,
                                        by         = interval,
                                        length.out = n + 1L),
                         rate = zeta))
    } else {
        h1l <- c(1.0, rep.int(0.0, n))
    }

    rp <- .C('rpbounds',
             h1l        = as.double(h1l),
             h1u        = as.double(c(0.0, h1l)),
             h2l        = as.double(h2l),
             h2u        = as.double(c(0.0, h2l)),
             q          = as.double(q),
             n          = as.integer(n),
             fl         = double(n + 3L),
             fu         = double(n + 3L),
             lowerbound = double(n + 3L),
             upperbound = double(n + 3L))

    ns <- seq_len(n)

    x <- seq.int(from       = 0.0,
                 by         = interval,
                 along.with = rp$lowerbound[ns])

    psi.lower <- stepfun(x     = x,
                         y     = c(rp$lowerbound[ns], 0.0),
                         right = TRUE)

    psi.upper <- stepfun(x     = x,
                         y     = c(1.0, rp$upperbound[ns]),
                         right = FALSE)

    if (use.splines) {
        psi <- splinefun(x = x,
                         y = 0.5 * (rp$lowerbound + rp$upperbound)[ns])
    } else {
        psi <- approxfun(x      = x,
                         y      = 0.5 * (rp$lowerbound + rp$upperbound)[ns],
                         method = 'linear',
                         rule   = 2L)
    }

    psi.x <- 0.5 * (rp$lowerbound + rp$upperbound)[ns]
    diff1 <- -diff(psi.x) / (interval * zeta * p)

    if (richardson) {
        diff3 <- -diff(psi.x, lag = 3L) / (3.0 * interval * zeta * p)
        diff1 <- rowSums(cbind(c(1.0, rep.int( 1.125, n - 3L)) * diff1[1L - n],    #  9 / 8
                               c(0.0, rep.int(-0.125, n - 3L)) * c(0.0, diff3)))   # -1 / 8
    }

    ns <- seq_len(n - 3L)

    x.1 <- c(0.0, x[ns] + 0.5 * interval)

    if (use.splines) {
        psi.1 <- splinefun(x = x.1,
                           y = c(1.0, diff1[ns]))
    } else {
        psi.1 <- approxfun(x      = x.1,
                           y      = c(1.0, diff1[ns]),
                           method = 'linear',
                           rule   = 2L)
    }

    return(structure(list(psi       = psi,
                          psi.1     = psi.1,
                          psi.2     = function(x) psi(x) - psi.1(x),
                          psi.lower = psi.lower,
                          psi.upper = psi.upper),
                     compmethod  = 'bounds',
                     riskproc    = process,
                     parameters  = list(interval   = interval,
                                        maxreserve = maxreserve),
                     diagnostics = list(rp = rp)))
}
