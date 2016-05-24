saddlepointRuinprob <- function(process, jensen = FALSE, normalize = TRUE) {
    stopifnot(is.logical(jensen), is.logical(normalize))

    p    <- process[['p']]
    zeta <- process[['zeta']]
    vx   <- process[['vx']]

    KL    <- process[['KL']]
    KL.d1 <- process[['KL.d1']]
    KL.d2 <- process[['KL.d2']]

    if (normalize) {
        corrconst <- integrate(
            f = function(v) {
                exp(KL(v) - v * KL.d1(v)) * sqrt(KL.d2(v)) * dnorm(0.0)
            },
            lower = -Inf,
            upper = adjcoef(process) - .Machine$double.eps^(2.0 / 3.0)
        )$value
    } else {
        corrconst <- 1.0
    }

    ## NOTE:  dnorm(0.0) == 1.0 / sqrt(2 * pi)

    if (jensen) { ## Jensen (1992)
        psi <- function(x) {
            res <- pnorm(vx(x)$z, lower.tail = FALSE)
            res[almost.equal(x, 0.0)] <- 1.0
            return(res)
        }

        psi.v <- function(v) {
            pnorm(process[['zv']](v), lower.tail = FALSE)
        }

        psi.1 <- function(x) {
            v <- vx(x)
            res <- pnorm(v$z - log(1.0 - v$v / (zeta * p * corrconst)) / v$r) - pnorm(v$z)
            res[almost.equal(x, 0.0)] <- 1.0
            return(res)
        }

        psi.1.v <- function(v) {
            z <- process[['zv']](v)
            r <- process[['rv']](v)
            pnorm(z - log(1.0 - v / (zeta * p * corrconst)) / r) - pnorm(z)
        }

        psi.2 <- function(x) {
            v <- vx(x)
            res <- pnorm(v$z - log(1.0 - v$v / (zeta * p * corrconst)) / v$r, lower.tail = FALSE)
            res[almost.equal(x, 0.0)] <- 0.0
            return(res)
        }

        psi.2.v <- function(v) {
            z <- process[['zv']](v)
            r <- process[['rv']](v)
            pnorm(z - log(1.0 - v / (zeta * p * corrconst)) / r, lower.tail = FALSE)
        }
    } else { # Lugannani and Rice (1980), Daniels (1954)
        psi <- function(x) {
            v <- vx(x)
            res <- pnorm(v$r, lower.tail = FALSE) - dnorm(v$r) * (1.0 / v$r - 1.0 / v$s)
            res[almost.equal(x, 0.0)] <- 1.0
            return(res)
        }

        psi.v <- function(v) {
            r <- process[['rv']](v)
            s <- process[['sv']](v)
            pnorm(r, lower.tail = FALSE) - dnorm(r) * (1.0 / r - 1.0 / s)
        }

        psi.1 <- function(x) {
            v <- vx(x)
            res <- exp(v$value) * dnorm(0.0) / (sqrt(v$hessian) * p * zeta * corrconst)
            res[almost.equal(x, 0.0)] <- 1.0
            return(res)
        }

        psi.1.v <- function(v) {
            exp(KL(v) - v * KL.d1(v)) / sqrt(KL.d2(v) * p * zeta * corrconst) * dnorm(0.0)
        }

        psi.2 <- function(x) {
            v <- vx(x)
            res <- pnorm(v$r, lower.tail = FALSE) - dnorm(v$r) * (1.0 / v$r - 1.0 / v$s * (1.0 - v$v / (p * zeta * corrconst)))
            res[almost.equal(x, 0.0)] <- 0.0
            return(res)
        }

        psi.2.v <- function(v) {
            r <- process[['rv']](v)
            s <- process[['sv']](v)
            pnorm(r, lower.tail = FALSE) - dnorm(r) * (1.0 / r - 1.0 / s * (1.0 - v / (p * zeta * corrconst)))
        }
    }

    return(structure(.Data       = list(psi     = psi,
                                        psi.1   = psi.1,
                                        psi.2   = psi.2),
                     compmethod  = 'saddlepoint',
                     riskproc    = process,
                     parameters  = list(jensen    = jensen,
                                        normalize = normalize),
                     diagnostics = list(corrconst = corrconst,
                                        psi.v     = psi.v,
                                        psi.1.v   = psi.1.v,
                                        psi.2.v   = psi.2.v)))
}
