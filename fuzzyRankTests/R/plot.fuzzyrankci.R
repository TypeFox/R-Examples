plot.fuzzyrankci <- function(x, y, add = FALSE, verticals = FALSE,
    col.hor = par("col"), col.vert = par("col"), lty.vert = 2,
    pch.vert = 19, full.ylim = TRUE, extra.xlim = 0.2, main, ...)
{
    if (missing(main))
        main <- "Fuzzy Confidence Interval"

    xx <- x$knots
    vv <- x$knot.values
    ii <- x$interval.values

    ##### make -Inf to Inf #####
    if (is.finite(min(x$knots))) {
        xx <- c(-Inf, xx)
        vv <- c(NA, vv)
        ii <- c(0, ii)
    }
    if (is.finite(max(x$knots))) {
        xx <- c(xx, Inf)
        vv <- c(vv, NA)
        ii <- c(ii, 0)
    }

    xxfin <- xx[is.finite(xx)]
    vvfin <- vv[is.finite(xx)]
    xxmin <- min(xxfin) - extra.xlim * diff(range(xxfin))
    xxmax <- max(xxfin) + extra.xlim * diff(range(xxfin))
    xxii <- xx
    xxii[1] <- xxmin
    xxii[length(xxii)] <- xxmax
    xxii <- (xxii[- length(xxii)] + xxii[- 1]) / 2
    xfoo <- c(xxfin, xxii)
    yfoo <- c(vvfin, ii)
    xfoo <- c(xxmin, xxmax, xfoo)
    yfoo <- c(min(yfoo), min(yfoo), yfoo)

    if (! add) {
        if (full.ylim)
            plot(xfoo, yfoo, type = "n", xlab = x$data.name,
                ylab = "", main = main, ylim = c(0, 1), ...)
        else
            plot(xfoo, yfoo, type = "n", xlab = x$data.name,
                ylab = "", main = main, ...)
    }

    fred <- par("usr")
    xx[1] <- fred[1]
    xx[length(xx)] <- fred[2]

    lenxx <- length(xx)
    x0 <- xx[- lenxx]
    x1 <- xx[- 1]
    segments(x0, ii, x1, ii, col = col.hor, ...)
    lenii <- length(ii)
    y0 <- ii[- lenii]
    y1 <- ii[- 1]
    if (verticals) {
        segments(xxfin, y0, xxfin, y1, col = col.vert, lty = lty.vert, ...)
    }
    points(xx, vv, pch = pch.vert, ...)
}
