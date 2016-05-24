
# C_xspline(x, y, s, open, repEnds, draw, col, border, ...)

# TODO: handle NA's in x and/or y (see polygon.R)
C_xspline <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:9)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    xx <- tx(x[[2]], par)
    yy <- ty(x[[3]], par)
    s <- x[[4]]
    open <- x[[5]]
    repEnds <- x[[6]]
    draw <- x[[7]]
    col <- FixupCol(x[[8]], NA, par$bg)
    border <- FixupCol(x[[9]], par$fg, par$bg)
    if (draw) {
        grid.xspline(xx, yy, default.units="native",
                     shape=s, open=open, repEnds=repEnds,
                     gp=gpar(col=border, fill=col, lwd=par$lwd, lty=par$lty,
                         lineend=par$lend, linemitre=par$lmitre,
                         linejoin=par$ljoin),
                     name=grobname("xspline"))
        result <- NULL
    }
    upViewport(depth)
}
