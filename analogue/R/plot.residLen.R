## plot method for "residLen"
`plot.residLen` <- function(x, probs = c(0.9, 0.95, 0.99),
                            ncol = 1, lcol = "red",
                            llty = "dashed",
                            xlab = NULL, ylab = NULL,
                            main = "Residual distances",
                            rug = TRUE,
                            ...) {
    ## compute densities
    d.train <- with(x, density(train, from = 0))
    d.passive <- with(x, density(passive, from = 0))
    y.lim <- range(0, d.train$y, d.passive$y)
    x.lim <- range(0, x$train)
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
        ylab <- "Density"
    }
    ## plotting
    op <- par(oma = c(4,0,5,0), las = 1, no.readonly = TRUE,
              mar = c(1.5,4,1,2) + 0.1)
    on.exit(par(op))
    layout(matrix(1:2, ncol = ncol))
    plot(d.train, type = "n", ann = FALSE, axes = FALSE,
         xlim = x.lim, ylim = y.lim)
    abline(h = 0, col = "lightgrey")
    abline(v = q.train, col = lcol, lty = llty)
    with(x, rug(train, side = 1))
    lines(d.train)
    axis(1)
    axis(2)
    axis(3, at = q.train,
         labels = paste(round(probs*100), "%"),
         las = 2)
    box()
    title(ylab = ylab)
    plot(d.passive, type = "n", ann = FALSE, axes = FALSE,
         xlim = x.lim, ylim = y.lim)
    abline(h = 0, col = "lightgrey")
    abline(v = q.train, col = lcol, lty = llty)
    with(x, rug(passive, side = 1))
    lines(d.passive)
    axis(1)
    axis(2)
    box()
    title(ylab = ylab)
    title(xlab = xlab, outer = TRUE, line = 2)
    title(main = main, outer = TRUE, line = 3)
    layout(1)
    invisible(list(train = d.train, passive = d.passive))
}
