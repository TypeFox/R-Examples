## MMD.R (2013-12-12)

##   Mismatch Distribution

## Copyright 2009-2013 Emmanuel Paradis, 2013 David Winter

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

MMD <- function(x, xlab = "Distance", main = "", rug = TRUE, legend = TRUE,
                lcol = c("blue", "red"), lty = c(1, 1), ...)
{
    d <- dist.dna(x, "N")
    hist(d, xlab = xlab, main = main, freq = FALSE, ...)
    lines(density(d, bw = 2), col = lcol[1], lty = lty[1])
    ## by David Winter:
    theta <- mean(d)
    upper <- ceiling(max(d))
    e <- sapply(0:upper, function(i) theta^i / (theta + 1)^(i + 1))
    lines(e, col = lcol[2], lty = lty[2])
    if (rug) rug(d)
    if (legend) {
        psr <- par("usr")
        xx <- psr[2]/2
        yy <- psr[4] * (0.5 + 0.5/par("plt")[4])
        legend(xx, yy, c("Empirical", "Stable expectation"),
               lty = lty, col = lcol, bg = "white", bty = "n",
               xjust = 0.5, yjust = 0.5, horiz = TRUE, xpd = TRUE)
        #legend("topleft", c("Empirical", "Stable expectation"),
        #       lty = 1, col = lcol, bg = "white", bty = "n")
    }
}
