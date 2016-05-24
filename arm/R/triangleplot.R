
triangleplot <- function (x, y = NULL, cutpts = NULL, details = TRUE, n.col.legend = 5,
    cex.col = 0.7, cex.var = 0.9, digits = 1, color = FALSE)
{
    if (!is.matrix(x))
        stop("x must be a matrix!")
    if (dim(x)[1] != dim(x)[2])
        stop("x must be a square matrix!")
    x.na <- x
    x.na[is.na(x.na)] <- -999
    z.plot <- x
    if (is.null(y)) {
        z.names <- dimnames(x)[[2]]
    }
    else {
        z.names <- y
    }
    for (i in 1:dim(z.plot)[1]) for (j in i:dim(z.plot)[2]) z.plot[i,
        j] <- NA
    layout(matrix(c(2, 1), 1, 2, byrow = FALSE), c(10.5, 1.5))
    layout(matrix(c(2, 1), 1, 2, byrow = FALSE), c(10.5, 1.5))
    if (is.null(cutpts)) {
        if (details) {
            neg.check <- abs(sum(z.plot[z.plot < 0], na.rm = T))
            if (neg.check > 0) {
                z.breaks <- sort(c(0, seq(min(z.plot, na.rm = T),
                  max(z.plot, na.rm = T), length = n.col.legend)))
            }
            else {
                z.breaks <- seq(min(z.plot, na.rm = T), max(z.plot,
                  na.rm = T), length = n.col.legend + 1)
            }
            for (i in 1:4) {
                n1 <- length(unique(round(z.breaks, digits = digits)))
                n2 <- length(z.breaks)
                ifelse((n1 != n2), digits <- digits + 1, digits <- digits)
            }
            if (digits > 3) {
                stop("Too many digits! Try to adjust n.col.legend to get better presentation!")
            }
        }
        else {
            postive.z <- na.exclude(unique(round(z.plot[z.plot >
                0], digits = digits)))
            neg.check <- abs(sum(z.plot[z.plot < 0], na.rm = T))
            ifelse(neg.check > 0, negative.z <- na.exclude(unique(round(z.plot[z.plot <
                0], digits = digits))), negative.z <- 0)
            max.z <- max(z.plot, na.rm = T)
            min.z <- min(z.plot, na.rm = T)
            z.breaks <- sort(unique(c(postive.z, negative.z)))
            n.breaks <- length(z.breaks)
            l.legend <- ceiling(n.col.legend/2)
            if (n.breaks > 8) {
                if (neg.check > 0) {
                  postive.z <- seq(0, max(postive.z), length = l.legend +
                    1)
                  negative.z <- seq(min(negative.z), 0, length = l.legend)
                  z.breaks <- sort(unique(c(postive.z, negative.z)))
                  n.breaks <- length(z.breaks)
                  z.breaks[1] <- min.z
                  z.breaks[n.breaks] <- max.z
                  n.col.legend <- length(z.breaks) - 1
                }
                else {
                  postive.z <- seq(0, max(postive.z), length = n.col.legend +
                    1)
                  z.breaks <- sort(unique(c(postive.z, negative.z)))
                  n.breaks <- length(z.breaks)
                  z.breaks[1] <- min.z
                  z.breaks[n.breaks] <- max.z
                  n.col.legend <- length(z.breaks) - 1
                }
            }
            else {
                if (neg.check > 0) {
                  z.breaks <- sort(c(0, seq(min(z.plot, na.rm = T),
                    max(z.plot, na.rm = T), length = n.col.legend)))
                }
                else {
                  z.breaks <- seq(min(z.plot, na.rm = T), max(z.plot,
                    na.rm = T), length = n.col.legend + 1)
                }
            }
        }
    }
    if (!is.null(cutpts)) {
        z.breaks = cutpts
        n.breaks <- length(z.breaks)
        n.col.legend <- length(z.breaks) - 1
    }
    if (color) {
        z.colors <- heat.colors(n.col.legend)[n.col.legend:1]
    }
    else {
        z.colors <- gray(n.col.legend:1/n.col.legend)
    }
    par(mar = c(0.5, 0.1, 2, 0.1), pty = "m")
    plot(c(0, 1), c(min(z.breaks), max(z.breaks)), type = "n",
        bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    for (i in 2:(length(z.breaks))) {
        rect(xleft = 0.5, ybottom = z.breaks[i - 1], xright = 1,
            ytop = z.breaks[i], col = z.colors[i - 1])
        text(x = 0.45, y = z.breaks[i - 1], labels = format(round(z.breaks[i -
            1], digits)), cex = cex.col, adj = 1, xpd = TRUE)
    }
    rect(xleft = 0.5, ybottom = z.breaks[length(z.breaks)], xright = 1,
        ytop = z.breaks[length(z.breaks)], col = z.colors[length(z.colors)])
    text(x = 0.45, y = z.breaks[length(z.breaks)], labels = format(round(z.breaks[length(z.breaks)],
        digits)), cex = cex.col, adj = 1, xpd = TRUE)
    par(mar = c(0.1, 0.1, 2, 0.1), pty = "m")
    image(x = 1:dim(z.plot)[1], y = 1:dim(z.plot)[2], z = z.plot,
        xaxt = "n", yaxt = "n", bty = "n", col = z.colors, breaks = z.breaks,
        xlim = c(-2, dim(z.plot)[1] + 0.5), ylim = c(-1, dim(z.plot)[2] +
            0.5), xlab = "", ylab = "")
    text(x = 1:dim(z.plot)[1], y = 1:dim(z.plot)[2], labels = z.names,
        cex = cex.var, adj = 1, xpd = TRUE)
    for (i in 1:dim(z.plot)[1]) {
        for (j in i:dim(z.plot)[2]) {
            if (x.na[i, j] == -999 & i != j)
                points(x = j, y = i, pch = "x", cex = 0.9)
        }
    }
}
