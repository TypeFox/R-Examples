### ===== actuar: An R Package for Actuarial Science =====
###
### Histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon

hist.grouped.data <-
    function(x, freq = NULL, probability = !freq,
             density = NULL, angle = 45, col = NULL, border = NULL,
             main = paste("Histogram of", xname), xlim = range(cj),
             ylim = NULL, xlab = xname, ylab, axes = TRUE,
             plot = TRUE, labels = FALSE, ...)
{
    ## Group boundaries are in the environment of 'x'
    cj <- eval(expression(cj), envir = environment(x))
    nj <- x[, 2]

    ## If any frequency is non finite, omit the group
    keep <- which(is.finite(nj))
    nj <- nj[keep]
    cj <- cj[c(1, keep + 1)]

    ## Some useful values
    n <- sum(nj)                        # total number of observations
    h <- diff(cj)                       # group widths
    dens <- nj/(n * h)                  # group "densities"

    ## Cannot plot histogram with infinite group
    if (any(is.infinite(cj)))
        stop("infinite group boundaries")

    ## The rest is taken from hist.default()
    xname <- paste(deparse(substitute(x), 500), collapse = "\n")
    equidist <- diff(range(h)) < 1e-07 * mean(h)
    if (is.null(freq))
    {
        freq <- if (!missing(probability))
            !as.logical(probability)
        else equidist
    }
    else if (!missing(probability) && any(probability == freq))
        stop("'probability' is an alias for '!freq', however they differ.")
    mids <- 0.5 * (cj[-1] + cj[-length(cj)])
    r <- structure(list(breaks = cj, counts = nj, intensities = dens,
                        density = dens, mids = mids, xname = xname,
                        equidist = equidist),
                   class = "histogram")
    if (plot)
    {
        plot(r, freq = freq, col = col, border = border, angle = angle,
             density = density, main = main, xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, axes = axes, labels = labels, ...)
        invisible(r)
    }
    else
        r
}
