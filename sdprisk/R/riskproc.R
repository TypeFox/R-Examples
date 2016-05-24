riskproc <- function(claims, premium, freq, variance) {
    stopifnot(is.numeric(premium), premium > 0, is.numeric(freq), freq > 0,
              is.numeric(variance), variance >= 0, "mu" %in% names(claims),
              is.claiminfo(claims))

    mu   <- mean(claims)
    p    <- 1.0 - freq * mu / premium
    q    <- freq * mu / premium
    zeta <- 2.0 * premium / variance

    if (is.hypoexp(claims)) {
        r <- attr(hypoexpRuinprob(list(p = p, q = q, zeta = zeta, claims = claims)), 'diagnostics')$r
        adjcoef <- min(Re(r[which(almost.equal(Im(r), 0.0) & Re(r) > 0.0)]))
    } else {
        mgfMaxAggregateLoss <-  function(x) {
            freq * claims$mgf(x) + 0.5 * variance * x^2.0 - premium * x - freq
        }
        # Find the largest root to get an upper bound of the interval
        adjcoef <- max(uniroot.all(f     = mgfMaxAggregateLoss,
                                   lower = .Machine$double.eps^0.25,
                                   upper = 1.0e+8))
        # Find the smallest non-zero root with high accuracy
        adjcoef <- min(uniroot.all(f     = mgfMaxAggregateLoss,
                                   lower = .Machine$double.eps^0.25,
                                   upper = adjcoef + 1.0,
                                   n     = 1.0e+5))
    }

    mgf    <- claims$mgf
    mgf.d1 <- claims$mgf.d1
    mgf.d2 <- claims$mgf.d2

    if (is.function(mgf)) {

        KL <- function(x) {
            is.na(x) <- (x > adjcoef)
            is.na(x) <- almost.equal(x, adjcoef)
            res <- p * x / (x - x^2.0 / zeta + q / mu * (1.0 - mgf(x)))
            res[which(almost.equal(x, 0.0))] <- 1.0
            return(log(res))
        }

        KL.osc <- function(x) {
            log1p(x * exp(KL(x)) / (p * zeta))
        }

        ML.osc <- function(x) {
            is.na(x) <- (x > adjcoef)
            is.na(x) <- almost.equal(x, adjcoef)
            res <- x^2.0 / (zeta * x - x^2.0 + zeta * q / mu * (1.0 - mgf(x))) + 1.0
            res[which(almost.equal(x, 0.0))] <- 1.0
            return(res)
        }

        if (is.function(mgf.d1)) {

            KL.d1 <- function(x) {
                is.na(x) <- (x > adjcoef)
                is.na(x) <- almost.equal(x, adjcoef)

                mgf.x    <- mgf(x)
                mgf.d1.x <- mgf.d1(x)

                res <- (x^2.0 * mu + q * zeta - q * zeta * mgf.x + x * mgf.d1.x * q
                        * zeta) / (x * zeta * mu - x^2.0 * mu + q * zeta - q * zeta *
                        mgf.x) / x

                mynan      <- (is.infinite(res) | is.nan(res) | almost.equal(x, 0.0))
                res[mynan] <- grad(KL, x[mynan])
                return(unlist(res))
            }

            if (is.function(mgf.d2)) {

                KL.d2 <- function(x) {
                    is.na(x) <- (x > adjcoef)
                    is.na(x) <- almost.equal(x, adjcoef)

                    mgf.x    <- mgf(x)
                    mgf.d1.x <- mgf.d1(x)
                    mgf.d2.x <- mgf.d2(x)

                    res <- (-(-2.0 * x * zeta^2.0 * mu * q * mgf.x + 2.0 * x^2.0 * zeta^2.0
                              * mu * mgf.d1.x * q + 4.0 * x^2.0 * mu * q * zeta * mgf.x - 4.0
                              * x^3.0 * mu * mgf.d1.x * q * zeta - 2.0 * q^2.0 * zeta^2.0
                              * mgf.x + q^2.0 * zeta^2.0 * mgf.x^2.0 - mgf.d2.x * q * zeta^2.0
                              * x^3.0 * mu + mgf.d2.x * q * zeta * x^4.0 * mu + x^2.0
                              * mgf.d2.x * q^2.0 * zeta^2.0 * mgf.x + 2.0 * x * zeta^2.0 * mu
                              * q - 4.0 * x^2.0 * mu * q * zeta - x^2.0 * mgf.d1.x^2.0 * q^2.0
                              * zeta^2.0 - x^2.0 * mgf.d2.x * q^2.0 * zeta^2.0 - x^4.0
                              * mu^2.0 + q^2.0 * zeta^2.0)
                            / ( - x * zeta * mu + x^2.0 * mu - q * zeta + q * zeta
                               * mgf.x)^2.0 / x^2.0)

                    mynan      <- (is.infinite(res) | is.nan(res) | almost.equal(x, 0.0))
                    res[mynan] <- grad(KL.d1, x[mynan])
                    return(unlist(res))
                }

            } else{
                KL.d2 <- NULL
            }
        } else {
            KL.d1 <- NULL
            KL.d2 <- NULL
        }
    } else {
        KL    <- NULL
        KL.d1 <- NULL
        KL.d2 <- NULL
    }

    if (!is.null(KL)) {

        vx <- function(x, vmin = -1.0e+64, vmax = adjcoef - .Machine$double.eps^0.65) {

            res <- Map(f        = optim,
                       arg      = x,
                       MoreArgs = list(par     = -1.0,
                                       fn      = function(y, arg) KL(y) - y * arg,
                                       gr      = function(y, arg) KL.d1(y) - arg,
                                       method  = 'BFGS',
                                       hessian = TRUE))

            v       <- vapply(res, `[[`, numeric(1L), 'par')
            value   <- vapply(res, `[[`, numeric(1L), 'value')
            hessian <- vapply(res, `[[`, numeric(1L), 'hessian')

            r <- sign(v) * sqrt(-2.0 * value)
            s <- v * sqrt(hessian)
            z <- r + log(s / r) / r

            return(list(v       = v,
                        r       = r,
                        s       = s,
                        z       = z,
                        value   = value,
                        hessian = hessian))
        }

        rv <- function(v) {
            sign(v) * sqrt(2.0 * (v * KL.d1(v) - KL(v)))
        }

        sv <- function(v) {
            v * sqrt(KL.d2(v))
        }

        zv <- function(v) {
            rv.v <- rv(v)
            rv.v + log(sv(v) / rv.v) / rv.v
        }

    }

    structure(.Data = list(premium   = premium,
                           freq      = freq,
                           variance  = variance,
                           diffusion = (variance != 0.0),
                           p         = p,
                           q         = q,
                           beta      = premium / (freq * mu) - 1.0,
                           zeta      = zeta,
                           KL        = KL,
                           KL.d1     = KL.d1,
                           KL.d2     = KL.d2,
                           KL.osc    = KL.osc,
                           ML.osc    = ML.osc,
                           vx        = vx,
                           rv        = rv,
                           sv        = sv,
                           zv        = zv,
                           claims    = claims,
                           adjcoef   = adjcoef),
              class = c('riskproc', 'list'))
}
