ruinprob <- function(x, param.list, compmethod = c('dg', 'exp'), flmethod = c('nonp', 'exp', 'lnorm', 'custom'), reserve, loading, fl = NA, interval = 0.5, implementation = c('R', 'C'), ...) {
    stopifnot(is.numeric(x))

    implementation <- match.arg(implementation)

    if (!missing(param.list) && is.list(param.list)) {
        try({
                compmethod <- param.list$compmethod
                flmethod   <- param.list$flmethod
                reserve    <- param.list$reserve
                loading    <- param.list$loading
                fl         <- param.list$fl
                interval   <- param.list$interval
            },
            silent = TRUE
        )
    }

    stopifnot(reserve  >= 0.0,
              loading  >= 0.0,
              interval >  0.0)

    if (is.array(x)) {
        apply(X          = x,
              MARGIN     = head(seq_along(dim(x)), -1L),
              #MARGIN     = seq_along(dim(x))[-1L],
              FUN        = ruinprob,
              compmethod = compmethod,
              flmethod   = flmethod,
              reserve    = reserve,
              loading    = loading,
              interval   = interval,
              fl         = fl,
              ...)
    } else {
        x <- as.vector(x)
        x <- x[is.finite(x)]

        compmethod <- match.arg(compmethod)
        flmethod   <- match.arg(flmethod)

        if (flmethod == 'custom') {
            stopifnot(is.function(fl))
        } else {
            if (is.function(fl)) {
                warning(paste("The option 'fl' is ignored for this method (flmethod = '", flmethod, "').", sep = ''))
            }
        }

        switch(compmethod,
            #exp = 1.0 / (1.0 + loading) * exp(-reserve * loading / (mean(x) * (1.0 + loading))),
            exp = dexp(reserve * loading / mean(x), 1.0 / (1.0 + loading)),
            dg  = {
                psi.0 <- 1.0 / (1.0 + loading)
                num   <- floor(reserve / interval) + 1L

                if(flmethod == "lnorm"){
                    lx      <- log(x)
                    mymu    <- mean(lx)
                    mysigma <- sd(lx)
                    myseq   <- seq.int(from = 0.0, by = interval, length.out = num)
                } else {
                    myseq   <- seq.int(from = 0.0, by = interval, length.out = num + 1L)
                }

                h.l <- switch(
                    flmethod,
                    nonp   = diff(sapply(X   = myseq,
                                         FUN = function(myarg){
                                             sum(pmin(myarg, x))
                                         })) / sum(x),
                    exp    = diff(pexp(myseq, 1.0 / mean(x))),
                    lnorm  = sapply(X   = myseq,
                                    FUN = function(x) {
                                        exp(-mymu - mysigma^2.0 / 2.0) *
                                        integrate(f          = plnorm,
                                                  lower      = x,
                                                  upper      = x + interval,
                                                  meanlog    = mymu,
                                                  sdlog      = mysigma,
                                                  lower.tail = FALSE)$value
                                    }),
                    custom = diff(sapply(myseq, fl, ...))
                )

                h.u <- c(0.0, h.l)
                #h.u <- c(0.0, head(h.l, -1L))

                if (reserve < interval) {
                    if (reserve == 0.0) {
                        return(psi.0)
                    } else {
                        #lower.limit <- 1.0 - f.l[1L]
                        #upper.limit <- psi.0
                        return(psi.0 * (1.0 + (1.0 - h.l[1L])/(1.0 - psi.0 * h.l[1L])) / 2.0)
                    }
                } else {
                    if (implementation == 'R') {

                        f.l <- vector('numeric', num)
                        f.u <- vector('numeric', num)

                        f.l[1L] <- (1.0 - psi.0) / (1.0 - psi.0 * h.l[1L])
                        f.u[1L] <- 1.0 - psi.0

                        fac.l <- psi.0 / (1.0 - psi.0 * h.l[1L])

                        for(i in tail(seq_len(num), -1L)){
                            #f.l[i] <- crossprod(h.l[2L:i], f.l[(i - 1L):1L]) * fac.l
                            #f.u[i] <- crossprod(h.u[2L:i], f.u[(i - 1L):1L]) * psi.0
                            f.l[i] <- crossprod(h.l[tail(seq_len(i), -1L)], f.l[rev(seq_len(i - 1L))]) * fac.l
                            f.u[i] <- crossprod(h.u[tail(seq_len(i), -1L)], f.u[rev(seq_len(i - 1L))]) * psi.0
                        }

                        lower.limit <- 1.0 - pmax(0.0, pmin(1.0, sum(head(f.l, -1L))))
                        upper.limit <- 1.0 - pmax(0.0, pmin(1.0, sum(f.u)))

                        return(mean(c(lower.limit, upper.limit), na.rm = TRUE))
                        
                    } else {
                        .C('rplimits',
                           h.l   = as.double(h.l),
                           #h.u  = as.double(h.u),
                           h.u   = as.double(c(0.0, head(h.l, -1L))),
                           psi.0 = as.double(psi.0),
                           num   = as.integer(num),
                           rp    = double(1L))$rp
                    }
                }
            }
        )
    }
}
