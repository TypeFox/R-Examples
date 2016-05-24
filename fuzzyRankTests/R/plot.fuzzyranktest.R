plot.fuzzyranktest <- function(x, y, type = c("pdf", "cdf"), add = FALSE,
    col.hor = par("col"), col.vert = par("col"), lty.vert = 2,
    extra.xlim = 0.2, main, ...)
{
    type <- match.arg(type)
    xx <- x$knots
    vv <- x$values

    xxmin <- max(0, min(xx) - extra.xlim * diff(range(xx)))
    xxmax <- min(1, max(xx) + extra.xlim * diff(range(xx)))
    xxfoo <- xx
    vvfoo <- vv
    if (xxmin < min(xx)) {
        xxfoo <- c(xxmin, xxfoo)
        vvfoo <- c(0, vvfoo)
    }
    if (xxmax > max(xx)) {
        xxfoo <- c(xxfoo, xxmax)
        vvfoo <- c(vvfoo, 1)
    }

    if (0 < min(xx)) {
        xx <- c(0, xx)
        vv <- c(0, vv)
    }
    if (1 > max(xx)) {
        xx <- c(xx, 1)
        vv <- c(vv, 1)
    }
    if (type == "pdf") {
        if (missing(main))
            main <- "PDF of Fuzzy P-value"
        xxbar <- rep(xxfoo, each = 2)
        vvbar <- c(0, rep(diff(vvfoo) / diff(xxfoo), each = 2), 0)
        if (! add)
            plot(xxbar, vvbar, type = "n", xlab = "significance level",
                ylab = "probability density", main = main, ...)
        ii <- diff(vv) / diff(xx)
        lenxx <- length(xx)
        x0 <- xx[- lenxx]
        x1 <- xx[- 1]
        segments(x0, ii, x1, ii, col = col.hor, ...)
        x0 <- xx[- c(1, lenxx)]
        y0 <- ii[- length(ii)]
        y1 <- ii[- 1]
        segments(x0, y0, x0, y1, lty = lty.vert, col = col.vert, ...)
    }
    if (type == "cdf") {
        if (missing(main))
            main <- "CDF of Fuzzy P-value"
        if (! add)
            plot(xxfoo, vvfoo, type = "n", xlab = "significance level",
                ylab = "probability", main = main, ...)
        lines(xx, vv, ...)
    }
}
