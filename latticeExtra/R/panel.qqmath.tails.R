##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##


panel.qqmath.tails <-
    function(x, f.value = NULL, distribution = qnorm,
    groups = NULL, ..., approx.n = 100, tails.n = 10)
{
    if (getRversion() >= "2.11.0") 
        .Deprecated(msg = paste("'panel.qqmath.tails' is deprecated.",
                    "Use 'panel.qqmath' from lattice 0.18-4 onwards."))
    x <- as.numeric(x)
    distribution <- if (is.function(distribution))
        distribution
    else if (is.character(distribution))
        get(distribution)
    else eval(distribution)
    nobs <- sum(!is.na(x))
    if (!is.null(groups)) {
        panel.superpose(x, y = NULL, f.value = f.value, distribution = distribution,
                        groups = groups, panel.groups = panel.qqmath.tails,
                        ..., approx.n = approx.n, tails.n = tails.n)
        return()
    }
    if (nobs == 0) return()
    pp <- ppoints(nobs)
    y <- sort(x)
    qq <- distribution(pp)
    if (nobs > approx.n + tails.n*2) {
        keep <- rep(FALSE, nobs)
        ## keep lowest and highest points
        keep[c(1:tails.n, nobs+1-(1:tails.n))] <- TRUE
        ## keep points spaced equally along distribution
        qqkeep <- seq(qq[tails.n+1], qq[nobs-tails.n], length=approx.n)
        keep[findInterval(qqkeep, qq)] <- TRUE
        qq <- qq[keep]
        y <- y[keep]
    }
    panel.xyplot(x = qq, y = y, ...)
}

