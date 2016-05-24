plotSignal <-
function (x, my.colors = c("black", "red", "blue"), ylab = "Peak Intensity",
    xlab = c("individuals", "Phenotype"), case.control = NULL,
    cex.leg = 0.6, dens.bw = "nrd0", dens.adjust = 1, n = 0,
    ...)
{
    foo <- function(...) {
        old.mfrow <- par("mfrow")
        old.mar <- par("mar")
        on.exit(par(mfrow = old.mfrow, mar = old.mar))
        mm <- matrix(c(2:1), nrow = 1, ncol = 2, byrow = TRUE)
        layout(mm, widths = c(1, 0.4))
        den <- density(x, bw = dens.bw, adjust = dens.adjust)
        par(mar = c(5.1, 0, 4.1, 2.1))
        plot(den$y, den$x, type = "l", axes = FALSE, xlab = "density",
            ylab = "", col = "red")
        polygon(den$y, den$x, col = "red1")
        ll <- par("usr")[3:4]
        par(mar = c(5.1, 4.1, 4.1, 2.1))
        if (!is.null(cutoffs)) {
            if (all(x <= max(cutoffs)))
                cat("WARNING! No data above maximum cutoff point\n")
            if (all(x >= min(cutoffs)))
                cat("WARNING! No data below minimum cutoff point\n")
            x.ord <- as.integer(cut(x, c(-Inf, cutoffs, Inf)))
        }
        else x.ord <- rep(1, length(x))
        if (is.null(case.control)) {
            plot(x, ylim = ll, yaxs = "i", xlab = xlab[1], type = "n", ylab = ylab, ...)
            points(x, col = my.colors[x.ord])
        }
        else {
            tt <- unique(case.control)
            if (length(tt) == 1) {
                stop("case.control must have 2 differents values at least")
            }
            if (length(tt) > 2) {
                plot(case.control, x, col = my.colors[x.ord],
                  ylim = ll, yaxs = "i", xlab = xlab[2], ylab = ylab,
                  ...)
            }
            if (length(tt) == 2) {
                plot(x, ylim = ll, yaxs = "i", xlab = xlab[1],
                  type = "n", ylab = ylab, ...)
                o <- case.control == tt[1]
                n <- sum(o)
                points(1:n, x[o], col = my.colors[x.ord[o]],
                  pch = 16)
                o <- case.control == tt[2]
                points((n + 1):(n + sum(o)), x[o], col = my.colors[x.ord[o]],
                  pch = 4)
                legend("bottomright", as.character(tt), pch = c(16,
                  4), title = "Case-control status", bty = "n",
                  horiz = TRUE, cex = cex.leg)
            }
        }
        if (!is.null(cutoffs)) {
            cutoffs <- sort(cutoffs)
            legend("bottomleft", legend = round(rev(cutoffs),
                4), bty = "n", lty = rev(1:length(cutoffs)),
                title = "Cut off points:", cex = cex.leg)
            abline(h = cutoffs, lty = 1:length(cutoffs))
        }
    }
    cutoffs = NULL
    foo()
    if (n > 0) {
        cat("Place cut off point(s) using locator...\n")
        cutoffs <- locator(n)$y
        foo()
    }
    invisible(cutoffs)
}
