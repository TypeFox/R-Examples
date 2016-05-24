
## ceeboo 2007

circleplot.dist <- function(x, cutoff = 0.5, col = 1, circle = FALSE, scale = 1.4) {
    if (!inherits(x, "dist"))
        stop("'x' not of class dist")

    x <- order.dist(x)          # seriation

    z  <- seq(-pi, pi, length.out = attr(x, "Size") + 1)
    x0 <- cos(z)
    y0 <- sin(z)

    r <- c(-1,1) * scale

    plot(x0, y0, type = "p", xlim = r, ylim = r, xlab = "", ylab = "", 
                 xaxt = "n", yaxt = "n", pty="s",)
    if (circle) {
        z <- seq(-pi, pi, 0.01)
        lines(cos(z), sin(z), lty = 2)
    }
    text(x0, y0, labels = dimnames(x), pos = sign(x0) + 3)

    w <- c(cut(c(x), seq(0, cutoff, length.out = 4)))

    k <- !is.na(w)

    if (any(k)) {
        i <- row.dist(x)[k]
        j <- col.dist(x)[k]

        segments(x0[i], y0[i], x0[j], y0[j], lwd = w[k], col = col)
    }
    invisible()
}

##
