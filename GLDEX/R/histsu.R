#  File src/library/graphics/R/hist.R
#  Part of the R package, http://www.R-project.org
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#  This function was written by R Core Team and modified by Steve Su

"histsu" <-
function (x, breaks = "Sturges", freq = NULL, probability = !freq, 
    include.lowest = TRUE, right = TRUE, density = NULL, angle = 45, 
    col = NULL, border = NULL, main = paste("Histogram of", xname), 
    xlim = range(breaks), ylim = NULL, xlab = xname, ylab, axes = TRUE, 
    plot = TRUE, labels = FALSE, nclass = NULL, ...) 
{
    if (!is.numeric(x)) 
        stop("'x' must be numeric")
    xname <- paste(deparse(substitute(x), 500), collapse = "\n")
    n <- length(x <- x[is.finite(x)])
    use.br <- !missing(breaks)
    if (use.br) {
        if (!missing(nclass)) 
            warning("'nclass' not used when 'breaks' is specified")
    }
    else if (!is.null(nclass) && length(nclass) == 1) 
        breaks <- nclass
    use.br <- use.br && (nB <- length(breaks)) > 1
    if (use.br) 
        breaks <- sort(breaks)
    else {
        if (!include.lowest) {
            include.lowest <- TRUE
            warning("'include.lowest' ignored as 'breaks' is not a vector")
        }
        if (is.character(breaks)) {
            breaks <- match.arg(tolower(breaks), c("sturges", 
                "fd", "freedman-diaconis", "scott"))
            breaks <- switch(breaks, sturges = nclass.Sturges(x), 
                "freedman-diaconis" = , fd = nclass.FD(x), scott = nclass.scott(x), 
                stop("unknown 'breaks' algorithm"))
        }
        else if (is.function(breaks)) {
            breaks <- breaks(x)
        }
        if (!is.numeric(breaks) || is.na(breaks) || breaks < 
            2) 
            stop("invalid number of 'breaks'")
        breaks <- pretty.su(range(x), nint = breaks)
        nB <- length(breaks)
        if (nB <= 1) 
            stop("hist.default: pretty.su.su() error, breaks=", format(breaks))
    }
    h <- diff(breaks)
    equidist <- !use.br || diff(range(h)) < 1e-07 * mean(h)
    if (!use.br && any(h <= 0)) 
        stop("'breaks' are not strictly increasing")
    if (is.null(freq)) {
        freq <- if (!missing(probability)) 
            !as.logical(probability)
        else equidist
    }
    else if (!missing(probability) && any(probability == freq)) 
        stop("'probability' is an alias for '!freq', however they differ.")
    diddle <- 1e-07 * stats::median(diff(breaks))
    fuzz <- if (right) 
        c(if (include.lowest) -diddle else diddle, rep.int(diddle, 
            length(breaks) - 1))
    else c(rep.int(-diddle, length(breaks) - 1), if (include.lowest) diddle else -diddle)
    fuzzybreaks <- breaks + fuzz
    h <- diff(fuzzybreaks)
    storage.mode(x) <- "double"
    storage.mode(fuzzybreaks) <- "double"
    bin <- cut(x, breaks, include.lowest = TRUE)
    counts <- tabulate(bin, length(levels(bin)))

    if (any(counts < 0)) 
        stop("negative 'counts'. Internal Error in C-code for \"bincount\"")
    if (sum(counts) < n) 
        stop("some 'x' not counted; maybe 'breaks' do not span range of 'x'")
    dens <- counts/(n * h)
    mids <- 0.5 * (breaks[-1] + breaks[-nB])
    r <- structure(list(breaks = breaks, counts = counts, intensities = dens, 
        density = dens, mids = mids, xname = xname, equidist = equidist), 
        class = "histogram")
    if (plot) {
        plot(r, freq = freq, col = col, border = border, angle = angle, 
            density = density, main = main, xlim = xlim, ylim = ylim, 
            xlab = xlab, ylab = ylab, axes = axes, labels = labels, 
            ...)
        invisible(r)
    }
    else r
}
