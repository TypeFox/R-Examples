## site.spectrum.R (2014-01-22)

##   Site Frequency Spectrum

## Copyright 2009-2014 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

site.spectrum <- function(x, folded = TRUE, outgroup = 1)
{
    if (is.list(x)) x <- as.matrix(x)
    n <- dim(x)[1]
    if (n == 1 || is.vector(x))
        stop("only one sequence in the data set")

    ss <- seg.sites(x)

    if (folded) {
        more.than.two.states <- 0L
        spectrum <- integer(floor(n/2))
        for (i in ss) {
            bfi <- base.freq(x[, i], freq = TRUE)
            bfi <- bfi[bfi > 0]
            if (length(bfi) > 2)
                more.than.two.states <- more.than.two.states + 1L
            else {
                j <- min(bfi)
                spectrum[j] <- spectrum[j] + 1L
            }
        }
        if (more.than.two.states)
            warning(paste(more.than.two.states,
                          "sites with more than two states were ignored"))
    } else { # unfolded spectrum
        anc <- x[outgroup, ss, drop = TRUE]
        outgroup.state <- anc %in% as.raw(c(24, 40, 72, 136))
        if (s <- sum(!outgroup.state)) {
            warning(paste(s, "sites with ambiguous state were ignored"))
            ss <- ss[outgroup.state]
        }
        spectrum <- apply(x[, ss], 2, function(y) sum(y[outgroup] != y[-outgroup]))
        spectrum <- tabulate(spectrum, nrow(x) - 1)
    }
    class(spectrum) <- "spectrum"
    attr(spectrum, "folded") <- folded
    spectrum
}

plot.spectrum <- function(x, col = "red", main = NULL, ...)
{
    if (is.null(main)) {
        main <- "Site Frequency Spectrum"
        main <- if (attr(x, "folded")) paste("Folded", main) else paste("Unfolded", main)
    }

    barplot(as.numeric(x), names.arg = 1:length(x), col = col, main = main, ...)
}
