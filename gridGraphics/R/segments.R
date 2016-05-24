
# C_segments(x0, y0, x1, y1, col=col, lty=lty, lwd=lwd, ...)

C_segments <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:8)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    x0 <- tx(x[[2]], par)
    y0 <- ty(x[[3]], par)
    x1 <- tx(x[[4]], par)
    y1 <- ty(x[[5]], par)
    col <- FixupCol(x[[6]], NA, par$bg)
    lty <- FixupLty(x[[7]], par$lty)
    lwd <- FixupLwd(x[[8]], par$lwd)
    grid.segments(x0, y0, x1, y1, default.units="native",
                 gp=gpar(col=col, lty=lty, lwd=lwd, 
                     lineend=par$lend, linemitre=par$lmitre,
                     linejoin=par$ljoin),
                 name=grobname("segments"))
    upViewport(depth)    
}
