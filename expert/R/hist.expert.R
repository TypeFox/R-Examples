### ===== expert =====
###
### Histogram of the aggregated expert distribution
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
###          Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

hist.expert <-
    function(x, freq = NULL, probability = !freq,
             density = NULL, angle = 45, col = NULL, border = NULL,
             main = paste("Histogram of" , xname), xlim = NULL,
             ylim = NULL, xlab = "x", ylab = expression(f(x)), axes = TRUE,
             plot = TRUE, labels = FALSE, ...)
{
    ## Keep name of object for the main title
    xname <- paste(deparse(substitute(x), 500), collapse = "\n")

    ## Bounds
    y <- x$probs
    x <- x$breaks

    ## If any probability is non finite, omit the group
    keep <- which(is.finite(y))
    y <- y[keep]
    x <- x[c(1, keep + 1)]

    ## Some useful values
    h <- diff(x)                        # group widths
    dens <- y/h                         # group "densities"

    ## Cannot plot histogram with infinite group
    if (any(is.infinite(x)))
        stop("infinite group boundaries")

    ## The rest is taken from hist.default()
    equidist <- diff(range(h)) < 1e-07 * mean(h)
    if (is.null(freq))
    {
        freq <- if (!missing(probability))
            !as.logical(probability)
        else equidist
    }
    else if (!missing(probability) && any(probability == freq))
        stop("'probability' is an alias for '!freq', however they differ.")
    mids <- 0.5 * (x[-1] + x[-length(x)])
    r <- structure(list(breaks = x, probs = y, intensities = dens,
                        density = dens, mids = mids, xname = xname,
                        equidist = equidist),
                   class = "histogram")
    if (plot)
    {
        plot(r, freq = freq, col = col, border = border, angle = angle,
             density = density, main = main, xlim = range(x), ylim = ylim,
             xlab = xlab, ylab = ylab, axes = axes, labels = labels, ...)
        invisible(r)
    }
    else
        r
}
