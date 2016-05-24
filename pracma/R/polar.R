##
##  p o l a r . R
##


polar <- function(t, r, type="l", 
                        col = "blue", grcol = "darkgrey", bxcol = "black",
                        main = "Polar Plot", add = FALSE, ...) {
    stopifnot(is.numeric(t), is.numeric(r))
    n <- length(t)
    if (!is.vector(t) || !is.vector(r) || length(r) != n)
        stop("Arguments 't', 'r' have to be vectors of equal length.")
    z <- cbind(t, r)

    # transform coordinates
    xy <- pol2cart(z)
    if (n == 1) dim(xy) <- c(1, 2)
    hy <- hypot(xy[, 1], xy[, 2])

    if (!add) {
        # grid circle coordinates
        drs <- pretty(c(0, hy), min.n = 3)
        phi <- deg2rad(seq(0, 360, by = 2))
        cx <- cos(phi)
        cy <- sin(phi)
        
        # plot grid circles
        mad <- max(abs(drs))
        plot(c(-mad, mad), c(-mad, mad), type = "n", asp = 1,
             axes = FALSE, main = main,
             xlab = "", ylab = "")      # ann = FALSE
        box(col = bxcol)
        for (dr in drs) {
            lines(dr*cx, dr*cy, col=grcol, lty=3)
        }

        s1 <- 0.5 * mad; s2 <- mad * sqrt(3)/2
        # grid lines
        lines(c(-mad, mad), c(0, 0), col=grcol, lty=3)
        lines(c(0, 0), c(-mad, mad), col=grcol, lty=3)
        lines(c(-s2, s2), c(-s1, s1), col=grcol, lty=3)
        lines(c(-s2, s2), c(s1, -s1), col=grcol, lty=3)
        lines(c(-s1, s1), c(-s2, s2), col=grcol, lty=3)
        lines(c(-s1, s1), c(s2, -s2), col=grcol, lty=3)

        # grid annotation
        for (dr in drs) {
            text(0, dr, as.character(dr), adj = c(0.5, 1), col = grcol, cex = 0.75)
            last <- drs[length(drs)]
            text(last, 0, "0", pos = 4, offset = 0.2, col = grcol, cex = 0.75, font=2)
            text(-last, 0, "180", pos = 2, offset = 0.2, col = grcol, cex = 0.75, font=2)
            text(0, last, "90", pos = 3, offset=0.2, col = grcol, cex = 0.75, font=2)
            text(0, -last, "270", pos = 1, offset = 0.3, col = grcol, cex = 0.75, font=2)
            text(s2, s1, "30", adj = c(0, 0), col = grcol, cex = 0.75)
            text(s1, s2, "60", adj = c(0, 0), col = grcol, cex = 0.75)
            text(-s1, s2, "120", adj = c(1, 0), col = grcol, cex = 0.75)
            text(-s2, s1, "150", adj = c(1, 0), col = grcol, cex = 0.75)
            text(-s2, -s1, "210", adj = c(1, 1), col = grcol, cex = 0.75)
            text(-s1, -s2, "240", adj = c(1, 1), col = grcol, cex = 0.75)
            text(s1, -s2, "300", adj = c(0, 1), col = grcol, cex = 0.75)
            text(s2, -s1, "330", adj = c(0, 1), col = grcol, cex = 0.75)
        }
    }

    # Plot the function (type can be 'l', 'p', or 'n')
    lines(xy[, 1], xy[, 2], type = type, col = col, ...)

    invisible()
}
