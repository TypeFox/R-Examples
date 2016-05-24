
# contour(x, y, z, levels, labels, labcex, drawlabels, method,
#         vfont, col, lty, lwd)

C_contour <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:13)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    xx <- tx(x[[2]], par)
    yy <- ty(x[[3]], par)
    zz <- x[[4]]
    levels <- x[[5]]
    drawlabels <- x[[8]]
    col <- FixupCol(x[[11]], NA, par$bg)
    col <- ifelse(is.na(col), par$col, col)
    lty <- FixupLty(x[[12]], par$lty)
    lty <- ifelse(is.na(lty), par$lty, lty)
    lwd <- FixupLwd(x[[13]], par$lwd)
    lwd <- ifelse(is.na(lwd), par$lwd, lwd)
    if (drawlabels) {
        warning("gridGraphics cannot emulate labels on contour lines")
    }
    clines <- contourLines(xx, yy, zz, levels=levels)
    if (length(clines)) {
        for (i in 1:length(clines)) {
            c <- clines[[i]]
            grid.lines(c$x, c$y, default.units="native",
                       gp=gpar(col=col, lty=lty, lwd=lwd,
                           lineend=par$lend, linemitre=par$lmitre,
                           linejoin=par$ljoin),
                       name=paste(grobname("contour"), i, sep="-"))
        }
    }
    upViewport(depth)    
}
