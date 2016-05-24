
## based on the stat_smooth() function in ggplot2 package.

panel.smoother <-
    function(x, y, form = y ~ x, method = "loess", ...,
             se = TRUE, level = 0.95, n = 100,
             col = plot.line$col, col.se = col,
             lty = plot.line$lty, lwd = plot.line$lwd,
             alpha = plot.line$alpha,
             alpha.se = 0.25, border = NA,
             ## ignored (do not pass to method()):
             subscripts, group.number, group.value,
             type, col.line, col.symbol, fill,
             pch, cex, font, fontface, fontfamily)
{
    plot.line <- trellis.par.get("plot.line")
    if (!missing(col.line))
        col <- col.line
    ## allow 'form' to be passed as the first argument
    missing.x <- missing(x)
    if (!missing.x && inherits(x, "formula")) {
        form <- x
        missing.x <- TRUE
    }
    ## use 'x' and 'y' if given
    ## otherwise try to find them in the formula environment
    if (missing.x)
        x <- environment(form)$x
    if (missing(y))
        y <- environment(form)$y
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 1) 
        return()
    x <- as.numeric(x)[ok]
    y <- as.numeric(y)[ok]
    mod <- do.call(method,
                   c(alist(form, data = list(x = x, y = y)),
                     list(...)))
    ## use the limits of the data, or panel limits, whichever is smaller
    lims <- current.panel.limits()
    xrange <- c(max(min(lims$x), min(x)), min(max(lims$x), max(x)))
    xseq <- seq(xrange[1], xrange[2], length = n)
    pred <- predict(mod, data.frame(x = xseq), se = se)
    if (se) {
        std <- qnorm(level/2 + 0.5)
        panel.polygon(x = c(xseq, rev(xseq)),
                      y = c(pred$fit - std * pred$se,
                      rev(pred$fit + std * pred$se)),
                      col = col.se, alpha = alpha.se, border = border)
        pred <- pred$fit
    }
    panel.lines(xseq, pred, col = col, alpha = alpha,
                lty = lty, lwd = lwd)
}
