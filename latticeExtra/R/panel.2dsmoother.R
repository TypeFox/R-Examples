
panel.2dsmoother <-
    function(x, y, z, subscripts = TRUE,
             form = z ~ x * y, method = "loess", ...,
             args = list(), n = 100)
{
    if (length(subscripts) == 0)
        return()
    ## allow 'form' to be passed as the first argument
    missing.x <- missing(x)
    if (!missing.x && inherits(x, "formula")) {
        form <- x
        missing.x <- TRUE
    }
    ## use 'x', 'y', 'z' if given
    ## otherwise try to find them in the formula environment
    if (missing.x)
        x <- environment(form)$x
    if (missing(y))
        y <- environment(form)$y
    if (missing(z))
        z <- environment(form)$z
    x <- x[subscripts]
    y <- y[subscripts]
    z <- z[subscripts]
    ok <- is.finite(x) & is.finite(y) & is.finite(z)
    if (sum(ok) < 1) 
        return()
    x <- as.numeric(x)[ok]
    y <- as.numeric(y)[ok]
    z <- as.numeric(z)[ok]
    mod <- do.call(method,
                   c(alist(form, data = list(x = x, y = y, z = z)),
                     args))
    ## use the limits of the data, or panel limits, whichever is smaller
    lims <- current.panel.limits()
    xrange <- c(max(min(lims$x), min(x)), min(max(lims$x), max(x)))
    yrange <- c(max(min(lims$y), min(y)), min(max(lims$y), max(y)))
    xseq <- seq(xrange[1], xrange[2], length = n)
    yseq <- seq(yrange[1], yrange[2], length = n)
    ## zseq <- seq(min(z), max(z), length = n)
    grid <- expand.grid(x = xseq, y = yseq)
    fit <- predict(mod, grid)
    panel.levelplot(x = grid$x, y = grid$y, z = fit, subscripts = TRUE,
                    ...)
}
