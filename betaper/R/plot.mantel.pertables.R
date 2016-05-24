`plot.mantel.pertables` <-
function (x, xlab = "Environmental distance", ylab = "Sorensen's similarity index", 
    pch = 19, ...) 
{
    layout(matrix(c(1:3, 3), 2, 2, byrow = TRUE), heights = c(1, 
        2))
    par(mar = c(2, 4, 4, 1), bty = "n")
    h1 <- hist(x$simulation$results[1, ], breaks = 25, main = "Mantel correlation values", 
        xlab = "", col = "lightblue", border = "darkred")
    par(new = TRUE)
    plot(density(x$simulation$results[1, ]), col = "red", lty = 3, 
        lwd = 2, main = "", xlab = "", ylab = "", axes = FALSE, 
        xlim = c(min(h1$breaks), max(h1$breaks)))
    abline(v = median(x$simulation$results[1, ]), lwd = 1.5, 
        col = "red")
    par(mar = c(2, 4, 4, 1), bty = "n")
    h2 <- hist(x$simulation$results[2, ], breaks = 25, main = "p-values", 
        xlab = "", col = "lightblue", border = "darkred")
    par(new = TRUE)
    plot(density(x$simulation$results[2, ]), col = "red", lty = 3, 
        lwd = 2, main = "", xlab = "", ylab = "", axes = FALSE, 
        xlim = c(min(h2$breaks), max(h2$breaks)))
    abline(v = median(x$simulation$results[2, ]), lwd = 1.5, 
        col = "red")
    par(mar = c(5, 4, 4, 2), bty = "l", bg = "transparent")
    plot(x$raw$vegdist ~ x$raw$env.dist, cex.lab = 1.2, cex.axis = 1.2, 
        col = "lightblue", main = "Correlation between distance matrices", 
        xlab = xlab, ylab = ylab, pch = pch, ...)
    par(new = TRUE)
    sapply(x$simulation$vegdist, function(y) lines(lowess(x$raw$env.dist, 
        y, f = 1), col = "lightblue", lwd = 1.5))
    lines(lowess(x$raw$env.dist, x$raw$vegdist, f = 1), col = "red", 
        lwd = 1.5, lty = 3)
}

