testPwrArgs <-
    function(d.primeA, d.prime0, sample.size, nsim, alpha)
    ## ncat, tau
{
    stopifnot(length(d.primeA) == 1L,
              is.numeric(d.primeA),
              d.primeA >= 0,
              is.finite(d.primeA))
    stopifnot(length(d.prime0) == 1L,
              is.numeric(d.prime0),
              d.prime0 >= 0,
              is.finite(d.prime0))
    stopifnot(length(sample.size) == 1L || length(sample.size) == 2L,
              all(is.numeric(sample.size)),
              all(sample.size > 0),
              all(is.finite(sample.size)))
    stopifnot(length(nsim) == 1L,
              is.numeric(nsim),
              nsim > 0,
              is.finite(nsim))
    stopifnot(length(alpha) == 1L,
              is.numeric(alpha),
              alpha >= 0,
              alpha <= 1)
    invisible()
}


dodPwr_internal <-
    function(d.primeA, d.prime0, sample.size, nsim = 1e3, alpha = 0.05,
             statistic=c("likelihood", "Wilcoxon", "Pearson", "Wald"),
             alternative = c("two.sided", "less", "greater"),
             tau=NULL, ...)
{
    ## Match arguements:
    stopifnot(!is.null(tau))
    stat <- match.arg(statistic)
    alt <- match.arg(alternative)

    ## Compute p-values:
    if(stat == "Wilcoxon") {
### NOTE: special case here to avoid fitting the DOD model.
        pvals <- vapply(seq_len(nsim), function(i) {
            data <- dodSim(tau=tau, d.prime=d.primeA,
                           sample.size=sample.size, method.tau="user.defined")
            nlev <- ncol(data) ## do we need to do something about the
            ## data here?
            wtest <- wilcox.test(rep.int(1:nlev, data[2, ]),
                                 rep.int(1:nlev, data[1, ]),
                                 alternative=alt, ...)
            unname(wtest$p.value)
        }, FUN.VALUE=numeric(1L))
    } else if(stat == "Wald") {
### FIXME: Do we need this special case here?
        ctrl <- dodControl(do.warn=FALSE)
        pvals <- vapply(seq_len(nsim), function(i) {
            data <- dodSim(tau=tau, d.prime=d.primeA,
                           sample.size=sample.size, method.tau="user.defined")
            fm <- dod_fit(data[1, ], data[2, ], control=ctrl)
            vcov <- fm$vcov
            std.err <- NA_real_
            if(!is.null(vcov) && all(is.finite(vcov)))
                std.err <- sqrt(diag(vcov)[ncol(vcov)])
            if(!is.na(std.err)) {
                stat.value <- (fm$d.prime - d.prime0) / std.err
                pval <- if(!is.finite(stat.value)) NA_real_ else
                normalPvalue(stat.value, alternative=alt)
            } else
                pval <- NA_real_
            pval
        }, FUN.VALUE=numeric(1L))
    } else { ## stat %in% c("likelihood", "Pearson")
        ctrl <- dodControl(get.vcov=FALSE, get.grad=FALSE, do.warn=FALSE)
        pvals <- vapply(seq_len(nsim), function(i) {
            data <- dodSim(tau=tau, d.prime=d.primeA,
                           sample.size=sample.size, method.tau="user.defined")
            fm <- try(dod(data[1, ], data[2, ], d.prime0=d.prime0,
                          alternative=alt, statistic=stat,
                          control=ctrl), silent=TRUE)
            if(inherits(fm, "try-error") || fm$optRes$convergence != 0)
                NA_real_ else fm$p.value
        }, FUN.VALUE=numeric(1L))
    }
    power <- mean(pvals < alpha, na.rm=TRUE)
    n <- sum(!is.na(pvals))
    se.power <-
        if(!is.finite(power) || power == 1 || power == 0) NA_real_
        else sqrt(power * (1 - power) / sqrt(n))
    attr(power, "se(power)") <- se.power
    attr(power, "n.used") <- n
    power
}

dodPwr <-
    function(d.primeA, d.prime0=0, ncat = 4, sample.size, nsim = 1e3,
             alpha = 0.05,
             method.tau=c("LR.max", "equi.prob", "se.min", "user.defined"),## max.X2?
             statistic=c("likelihood", "Wilcoxon", "Pearson", "Wald"),
             alternative = c("difference", "similarity", "two.sided",
             "less", "greater"),
             tau=NULL, ...)
{
    ## Match arguments:
    method <- match.arg(method.tau)
    stat <- match.arg(statistic)
    alt <- match.arg(alternative)

    ## Basic tests of arguments:
    testPwrArgs(d.primeA, d.prime0, sample.size, nsim, alpha)
    if(abs(nsim - round(nsim)) > 1e-6)
        warning("non-integer 'nsim' rounded to ", round(nsim))
    nsim <- round(nsim)
    if(max(abs(sample.size - round(sample.size))) > 1e-6)
        warning("non-integer 'sample.size' rounded to ", round(sample.size))
    size <- round(sample.size)

    ## Test that d.prime0, d.primeA conforms with 'alternative'
    if(alt %in% c("difference", "greater") && d.primeA < d.prime0)
        stop("Need d.primeA >= d.prime0 when alternative is '", alt, "'")
    if(alt %in% c("similarity", "less") && d.primeA > d.prime0)
        stop("Need d.primeA <= d.prime0 when alternative is '", alt, "'")
    if(alt == "difference") alt <- "greater"
    if(alt == "similarity") alt <- "less"

    if(!isTRUE(all.equal(d.prime0, 0)) && stat == "Wilcoxon") {
        stop("Wilcoxon statistic only available with d.prime0 = 0")
    }

    ## Match and test tau-related args:
    tau <- get_tau(d.prime=d.primeA, tau=tau, ncat=ncat,
                   method.tau=method, d.prime0=d.prime0)

    ## Return:
    dodPwr_internal(d.primeA=d.primeA, d.prime0=d.prime0,
                    sample.size=size, nsim=nsim, alpha=alpha,
                    statistic=stat, alternative=alt, tau=tau, ...)
}


## dodPvals_LR <-
##     function(d.primeA, d.prime0, tau, nsim, size, alternative)
## {
##     pvals <- vapply(seq_len(nsim), function(i) {
##         data <- dodSim(tau=tau, d.prime=d.primeA, sample.size=size)
##         fm1 <- try(dod_fit(data[1, ], data[2, ], get.vcov=FALSE,
##                            get.grad=FALSE), silent=TRUE)
##         fm0 <- try(dod_fit(data[1, ], data[2, ], d.prime=d.prime0,
##                            get.vcov=FALSE, get.grad=FALSE), silent=TRUE)
##         pval <- NA_real_
##         if(!inherits(fm1, "try-error") && !inherits(fm0, "try-error")) {
##             conv1 <- fm1$optRes$convergence
##             conv1 <- if(!is.null(conv1)) conv1 == 0 else TRUE
##             conv0 <- fm0$optRes$convergence
##             conv0 <- if(!is.null(conv0)) conv0 == 0 else TRUE
##             if(conv1 && conv0) {
##                 LR <- 2 * (fm1$logLik - fm0$logLik) ## LR shold be positive
##                 LR[LR < 0 && abs(LR) < 1e-4] <- 0
##                 if(LR >= 0) {
##                     lroot <- sign(d.primeA - d.prime0) * sqrt(LR)
##                     pval <- unname(normalPvalue(lroot, alternative=alt))
##                 }
##             }
##         }
##         pval
##     }, FUN.VALUE=numeric(1L))
##     pvals
## }

