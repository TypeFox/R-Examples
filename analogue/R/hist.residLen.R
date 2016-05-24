`hist.residLen` <- function(x, breaks = "Sturges",
                            freq = TRUE,
                            probs = c(0.9, 0.95, 0.99),
                            ncol = 1, lcol = "red",
                            llty = "dashed",
                            xlab = NULL, ylab = NULL,
                            main = "Residual distances",
                            rug = TRUE,
                            ...) {
    ## quantiles of resid distance
    q.train <- with(x, quantile(train, probs = probs))
    ## labels
    if(is.null(xlab)) {
        if(attr(x, "method") == "cca") {
            xlab <- expression(paste(Squared ~ chi^2 ~
                residual ~ distance))
        } else {
            xlab <- "Squared Euclidean residual distance"
        }
    }
    if(is.null(ylab)) {
        if(freq) {
            ylab <- "Frequency"
        } else {
            ylab <- "Density"
        }
    }
    ## plotting
    op <- par(oma = c(4,0,5,0), las = 1, no.readonly = TRUE,
              mar = c(1.5,4,1,2) + 0.1)
    on.exit(par(op))
    layout(matrix(1:2, ncol = ncol))
    h.train <- hist(x$train, breaks = breaks, plot = FALSE, ...)
    h.passive <- hist(x$passive, breaks = h.train$breaks,
                      plot = FALSE, ...)
    y <- if (freq)
        c(h.train$counts, h.passive$counts)
    else {
        y <- c(h.train$density, h.passive$density)
        if (is.null(y))
            c(h.train$intensities, h.passive$intensities)
        else y
    }
    y.lim <- range(0, y)
    plot(h.train, main = "", xlim = range(h.train$breaks),
         ylab = ylab, freq = freq, ylim = y.lim)
    abline(v = q.train, col = lcol, lty = llty)
    with(x, rug(train, side = 1))
    axis(1)
    axis(2)
    axis(3, at = q.train,
         labels = paste(round(probs*100), "%"),
         las = 2)
    box()
    title(ylab = ylab)
    plot(h.passive, main = "", xlim = range(h.train$breaks),
         ylab = ylab, freq = freq, ylim = y.lim)
    abline(v = q.train, col = lcol, lty = llty)
    with(x, rug(passive, side = 1))
    axis(1)
    axis(2)
    box()
    title(ylab = ylab)
    title(xlab = xlab, outer = TRUE, line = 2)
    title(main = main, outer = TRUE, line = 3)
    layout(1)
    invisible(list(train = h.train, passive = h.passive))
}
