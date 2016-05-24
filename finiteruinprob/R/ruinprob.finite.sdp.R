ruinprob.finite.sdp <- function(mgf, mgf.d1, mgf.d2, premium, freq, variance, endpoint, verbose = FALSE) {

    ## If the right boundary point of the domain of k is not supplied, take
    ## a wild guess (assume some large value such that no parts of the actual
    ## domain should get missed).  However, this could cause wrong results,
    ## because optimization and root-finding algorithms could stumble into areas
    ## outside the domain of certain functions.  Thus, please avoid this
    ## situation.
    if (missing(endpoint) || is.null(endpoint) || !is.finite(endpoint)) {
        endpoint <- 1.0e+6
        warning(sQuote('endpoint'), ' missing or invalid, using ',
                endpoint, ' as a substitute.')
    }

    ## If the function passed as mgf argument has already a 'difforder'
    ## argument, assume the user did not do this by mistake!
    if ('difforder' %in% names(formals(mgf))) {
        MX <- mgf
    } else {
        ## Make the not-all-too-bizarre assumption that the second-order
        ## derivative is only available if the first-order one is too.
        whichmissing <- c(missing(mgf.d1), missing(mgf.d2))
        if (any(whichmissing)) {
            if (whichmissing[1L]) {
                mgf.d1 <- function(x) {
                    drop(diag(numDeriv::grad(mgf, x)))
                }
                mgf.d2 <- function(x) {
                    drop(diag(numDeriv::hessian(mgf, x)))
                }
            } else if (whichmissing[2L]) {
                mgf.d2 <- function(x) {
                    drop(diag(numDeriv::grad(mgf.d1, x)))
                }
            }
        }

        MX <- function(x, difforder = 0L) {
            switch(match(difforder, 0L:2L, 4L),
                   mgf(x),
                   mgf.d1(x),
                   mgf.d2(x),
                   rep.int(NA_real_, length(x)))
        }
    }


    ## Additional parameters.
    .mu   <- MX(0.0, 1L)
    .p    <- 1.0 - freq * .mu / premium
    .q    <- freq * .mu / premium
    .beta <- premium / (freq * .mu) - 1.0
    .zeta <- 2.0 * premium / variance
    eps   <- 0.0001


    ## CGF of the aggregate loss process at time 1 (and its first two
    ## derivatives).
    k <- function(x, difforder = 0L) {
        is.na(x) <- vapply(x > endpoint, isTRUE, logical(1L))

        switch(match(difforder, 0L:2L, 4L),
               ## 1
               0.5 * variance^2.0 * x^2.0 - premium * x + freq * (MX(x, 0L) - 1.0),
               ## 2
               variance^2.0 * x - premium + freq * MX(x, 1L),
               ## 3
               variance^2.0 + freq * MX(x, 2L),
               ## 4
               rep.int(NA_real_, length(x)))
    }


    ## This gives the minimum and the value at the minimum of the function k.
#     beta.zero  <- nlm(f = function(x) {
#                           structure(k(x, 0L),
#                                     gradient = k(x, 1L),
#                                     hessian  = k(x, 2L))
#                       },
#                       p = 0.0)
# 
#     alpha.zero <- -beta.zero$minimum
#     beta.zero  <- beta.zero$estimate

    beta.zero  <- optim(par    = 0.0,
                        fn     = function(x) k(x, 0L),
                        gr     = function(x) k(x, 1L),
                        method = 'BFGS')[c('par', 'value')]

    alpha.zero <- -beta.zero$value
    beta.zero  <-  beta.zero$par


    ## Calculate the adjustment coefficient as it will be used below.
    rtmp     <- uniroot(f     = k,
                        lower = beta.zero,
                        upper = endpoint - eps^2.0)$root

    .adjcoef <- uniroot(f     = k,
                        lower = 0.9 * rtmp,
                        upper = min(1.1 * rtmp, endpoint - eps^2.0),
                        tol   = .Machine$double.eps^0.5)$root

    ## The CGF of the maximal aggregate loss
    KL <- function(x, difforder = 0L) {
        difforder <- match(difforder, 0L:2L, 4L)

        calc <- function(x) {
            is.na(x) <- (x > .adjcoef - eps)

            mgf.x <- MX(x, 0L)
            if (difforder >= 2L) {
                mgf.d1.x <- MX(x, 1L)
            }
            if (difforder >= 3L) {
                mgf.d2.x <- MX(x, 2L)
            }

            switch(difforder,
                   {# 1
                       res <- .p * x / (x - x^2.0 / .zeta + .q / .mu * (1.0 - mgf.x))
                       is.na(res) <- is.nan(res)
                       is.na(res) <- (res < 0)
                       log(res)
                   },
                   {# 2
                       (x^2.0 * .mu + .q * .zeta - .q * .zeta * mgf.x + x * mgf.d1.x * .q
                        * .zeta) / (x * .zeta * .mu - x^2.0 * .mu + .q * .zeta - .q * .zeta
                        * mgf.x) / x
                   },
                   {# 3
                       -(-2.0 * x * .zeta^2.0 * .mu * .q * mgf.x + 2.0 * x^2.0 * .zeta^2.0
                         * .mu * mgf.d1.x * .q + 4.0 * x^2.0 * .mu * .q * .zeta * mgf.x
                         - 4.0 * x^3.0 * .mu * mgf.d1.x * .q * .zeta - 2.0 * .q^2.0
                         * .zeta^2.0 * mgf.x + .q^2.0 * .zeta^2.0 * mgf.x^2.0 - mgf.d2.x
                         * .q * .zeta^2.0 * x^3.0 * .mu + mgf.d2.x * .q * .zeta * x^4.0
                         * .mu + x^2.0 * mgf.d2.x * .q^2.0 * .zeta^2.0 * mgf.x + 2.0
                         * x * .zeta^2.0 * .mu * .q - 4.0 * x^2.0 * .mu * .q * .zeta
                         - x^2.0 * mgf.d1.x^2.0 * .q^2.0 * .zeta^2.0 - x^2.0 * mgf.d2.x
                         * .q^2.0 * .zeta^2.0 - x^4.0 * .mu^2.0 + .q^2.0 * .zeta^2.0) /
                       (-x * .zeta * .mu + x^2.0 * .mu - .q * .zeta + .q * .zeta
                        * mgf.x)^2.0 / x^2.0
                   },
                   {# 4
                       NA_real_
                   })
        }

        exactval  <- c(0.0, rep.int(NA_real_, 3L))[difforder]
        .eps      <- eps * c(1.0, 2.0, 10.0, 10.0)[difforder]
        badx      <- vapply(abs(x) < .eps, isTRUE, logical(1L))
        val       <- calc(x)
        val[badx] <- approx(x    = c(-.eps, .eps, 0.0),
                            y    = c(calc(c(-.eps, .eps)), exactval),
                            xout = x[badx])$y

        return(val)
    }


    ## The (improper) two-dimensional common CGF of the time to ruin (first
    ## component) and the initial reserve (second component), shifted by 1 unit
    ## in the second component and transformed in the first component.
    K.tilde <- function(x) {
        ## Are both arguments zero?
        # check <- c(isTRUE(all.equal(x[1L], 0.0)), isTRUE(all.equal(x[2L], 0.0)))
        check <- vapply(Map(all.equal, x, c(0.0, 0.0)), isTRUE, logical(1L))
        if (all(check)) {
            res <- -0.5 * k(0.0, 2L) / k(0.0, 1L)
        } else {
            ## Are both arguments the same?
            # if (isTRUE(all.equal(x[1L], x[2L]))) {
            if (isTRUE(Reduce(all.equal, x))) {
                res <- (k(x[1L], 0L) / (x[1L] * k(x[1L], 1L)) - 1.0) / x[1L]
            } else {
                kval <- k(x, 0L) * c(-1.0, 1.0)
                nume <- -kval / x
                ## Is one of the arguments zero?
                if (Reduce(xor, check)) {
                    nume[check] <- k(x[check], 1L) * c(1.0, -1.0)[check]
                }
                res <- sum(nume) / sum(kval)
            }
        }
        is.na(res) <- (res <= 0.0) || !is.finite(res)
        return(log(res))
    }

    K.tilde.opt <- function(x, target) {
        K.tilde(x) - drop(crossprod(c(-k(x[1L], 0L), x[2L]), target))
    }


    ## The saddlepoint approximation to the probability of ruin within a finite
    ## time horizon
    ##
    ## TODO:  It should be possible to use this function for t = Inf
    psi <- function(x, t) {
        target <- c(t, x)
        res <- list(optim(par     = c(-0.5, -0.5),
                          fn      = K.tilde.opt,
                          gr      = NULL,
                          target  = target,
                          method  = 'BFGS',
                          hessian = TRUE),
                    optim(par     = -0.5,
                          fn      = function(x) K.tilde.opt(c(0, x), target),
                          gr      = NULL,
                          method  = 'BFGS',
                          hessian = TRUE),
                    optim(par     = 0.0,
                          fn      = function(x) KL(x, 0L) - x * target[2L],
                          gr      = function(x) KL(x, 1L) - target[2L],
                          method  = 'BFGS',
                          hessian = TRUE))

        names(res) <- c('bivariate', 'condition', 'infinite')
        res        <- lapply(res, `[`,  c('par', 'value', 'hessian'))

        res[[c('bivariate', 'hessian')]]   <- det(res[[c('bivariate', 'hessian')]] / ((-k(res[[c('bivariate', 'par')]][1L], 1L))^matrix(c(2.0, 1.0, 1.0, 0.0), 2L, 2L)))
        res[[c('bivariate', 'par')]][[1L]] <- -k(res[[c('bivariate', 'par')]][1L], 0L)
        res[[c('condition', 'hessian')]]   <- drop(res[[c('condition', 'hessian')]])
        res[[c('infinite', 'hessian')]]    <- drop(res[[c('infinite', 'hessian')]])

        s.fin <- sqrt(res[[c('bivariate', 'hessian')]] / res[[c('condition', 'hessian')]]) * res[[c('bivariate', 'par')]][1L]
        r.fin <- sign(res[[c('bivariate', 'par')]][1L]) * sqrt(2.0 * (res[[c('condition', 'value')]] - res[[c('bivariate', 'value')]]))
        z.fin <- r.fin + log(s.fin / r.fin) / r.fin

        psi.fin      <- pnorm(r.fin) + (1.0 / r.fin - 1.0 / s.fin) * dnorm(r.fin)
        psi.fin.star <- pnorm(z.fin)

        s.inf <- res[[c('infinite', 'par')]] * sqrt(res[[c('infinite', 'hessian')]])
        r.inf <- sign(res[[c('infinite', 'par')]]) * sqrt(-2.0 * res[[c('infinite', 'value')]])
        z.inf <- r.inf + log(s.inf / r.inf) / r.inf

        psi.inf      <- pnorm(r.inf, lower.tail = FALSE) - (1.0 / r.inf - 1.0 / s.inf) * dnorm(r.inf)
        psi.inf.star <- pnorm(z.inf, lower.tail = FALSE)

        psi      <- psi.inf      * psi.fin
        psi.star <- psi.inf.star * psi.fin.star

        if (verbose) {
            return(structure(.Data = list(psi          = psi,
                                          psi.star     = psi.star),
                             inf   = list(psi.inf      = psi.inf,
                                          psi.inf.star = psi.inf.star),
                             fin   = list(psi.fin      = psi.fin,
                                          psi.fin.star = psi.fin.star),
                             diag  = list(sdp   = res,
                                          s.inf = s.inf,
                                          r.inf = r.inf,
                                          z.inf = z.inf,
                                          s.fin = s.fin,
                                          r.fin = r.fin,
                                          z.fin = z.fin)))
        } else {
            return(list(psi      = psi,
                        psi.star = psi.star))
        }
    }

    if (verbose) {
        return(structure(.Data = psi,
                         diag  = list(k           = k,    ## The CGF of the loss process at time 1
                                      KL          = KL,   ## The CGF of the maximal aggregate loss
                                      K.tilde     = K.tilde,
                                      K.tilde.opt = K.tilde.opt,
                                      premium     = premium,
                                      freq        = freq,
                                      variance    = variance,
                                      endpoint    = endpoint,
                                      alpha.zero  = alpha.zero,
                                      beta.zero   = beta.zero,
                                      adjcoef     = .adjcoef,
                                      eps         = eps)))
    } else {
        return(psi)
    }
}

### vim: tw=78 spell spelllang=en_gb
