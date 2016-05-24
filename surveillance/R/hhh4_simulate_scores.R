################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Compute scores based on simulations from fitted hhh4() models
###
### Copyright (C) 2013-2015 Sebastian Meyer
### $Revision: 1476 $
### $Date: 2015-09-15 00:08:30 +0200 (Die, 15. Sep 2015) $
################################################################################


## logarithmic score
## CAVE: will be infinite if none of "sims" yields "x"
logs_sims <- function (sims, x)
    .logs(px = mean(sims == x))

## Dawid-Sebastiani score
## CAVE: undefined if all simulations have the same value (i.e., no variance)
dss_sims <- function (sims, x)
{
    if ((varsims <- var(sims)) == 0) { # FIXME: What to do in that case?
        warning("DSS undefined for zero variance of prediction: all(sims==",
                sims[1L], "), x=", x)
        NA_real_ # if (x==sims[1L]) -Inf else Inf
    } else {
        .dss(meanP = mean(sims), varP = varsims, x = x)
    }
}

## ranked probability score
rps_sims <- function (sims, x)
{
    .rps(P = ecdf(sims), x = x, kmax = ceiling(mean(sims) + 40*sd(sims)))

    ## Two alternatives via the expectation-based definition of the RPS:
    ## method = "means": equivalent to ecdf approach but slower
    ## method = "means.MC": faster than ecdf but with approximation error
    ## simdiffs <- switch(method,
    ##                    "means.MC" = diff(sims),
    ##                    "means" = outer(sims, sims, "-"))
    ## mean(abs(sims - x)) - mean(abs(simdiffs)) / 2
}


## scores-method for simulations from a hhh4 fit
scores.hhh4sims <- function (x, which = "rps", units = NULL, ..., drop = TRUE)
{
    observed <- observed(attr(x, "stsObserved"))
    scoreFUNs <- mget(paste0(which, "_sims"),
                      envir = getNamespace("surveillance"), inherits = FALSE)
    names(scoreFUNs) <- which
    if (!is.null(units)) {
        observed <- observed[, units, drop = FALSE]
        x <- x[, units, , drop = FALSE]
    }
    
    counts <- array(c(observed, x), dim = dim(x) + c(0L, 0L, 1L))
    res <- lapply(X = scoreFUNs, FUN = function (scoreFUN)
        apply(counts, 1:2, function (y) scoreFUN(y[-1L], y[1L])))
    res <- simplify2array(res, higher = TRUE)
    if (drop) drop(res) else res
}

## scores-method for simulations from a bunch of hhh4 fits
scores.hhh4simslist <- function (x, ...)
    lapply(X = x, FUN = scores.hhh4sims, ...)
