"image_plot" <- function (..., add = FALSE, nlevel = 64, legend.shrink = 0.9, 
    legend.width = 0.05, graphics.reset = FALSE, horizontal = FALSE, 
    offset = 2 * legend.width, bigplot = NULL, smallplot = NULL, 
    legend.only = FALSE, col = topo.colors(nlevel)) {

  # this is modified slightly from the fields function image.plot() and uses the fields auxiliary functions image.plot.info() and image.plot.plt()
    old.par <- par(no.readonly = TRUE)
    info <- image_plot_info(...)
    if (add) 
        big.plot <- old.par$plt
    if (legend.only) 
        graphics.reset <- TRUE
    temp <- image_plot_plt(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, horizontal = horizontal, 
        offset = offset, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        image(..., add = add, col = col)
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    temp <- list(...)
    iy <- seq(info$zlim[1], info$zlim[2], , nlevel)
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    ix <- 1
    if (!horizontal) {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
        axis(4, mgp = c(3, 1, 0), las = 2)
        box()
    }
    else {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        image(iy, ix, t(iz), yaxt = "n", xlab = "", ylab = "", 
            col = col)
        box()
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}
