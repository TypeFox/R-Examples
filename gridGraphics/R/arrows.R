
# arrows(x0, y0, x1, y1, length, angle, code, col, lty, lwd, ...) 

C_arrows <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:11)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    x0 <- tx(x[[2]], par)
    y0 <- ty(x[[3]], par)
    x1 <- tx(x[[4]], par)
    y1 <- ty(x[[5]], par)
    length <- x[[6]]
    angle <- x[[7]]
    code <- x[[8]]
    col <- FixupCol(x[[9]], NA, par$bg)
    lty <- FixupLty(x[[10]], par$lty)
    lwd <- FixupLwd(x[[11]], par$lwd)
    grid.segments(x0, y0, x1, y1, default.units="native",
                  gp=gpar(col=col, lty=lty, lwd=lwd,
                      lineend=par$lend, linemitre=par$lmitre,
                      linejoin=par$ljoin),
                  arrow=arrow(angle=angle, length=unit(length, "in"),
                      ends=switch(code, "first", "last", "both")),
                  name=grobname("arrows"))
    upViewport(depth)
}
