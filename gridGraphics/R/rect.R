
# C_rect(xleft, ybottom, xright, ytop, col, border, lty, lwd, ...)

C_rect <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:9)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    xleft <- tx(x[[2]], par)
    ybottom <- ty(x[[3]], par)
    xright <- tx(x[[4]], par)
    ytop <- ty(x[[5]], par)
    col <- FixupCol(x[[6]], NA, par$bg)
    border <- FixupCol(x[[7]], par$fg, par$bg)
    lty <- FixupLty(x[[8]], par$lty)
    lty <- ifelse(is.na(lty), par$lty, lty)
    lwd <- FixupLwd(x[[9]], par$lwd)
    lwd <- ifelse(is.na(lwd), par$lwd, lwd)
    grid.rect(xleft, ybottom, xright - xleft, ytop - ybottom,
              default.units="native", just=c("left", "bottom"),
              gp=gpar(col=border, fill=col, lty=lty, lwd=lwd,
                  lineend=par$lend, linemitre=par$lmitre,
                  linejoin=par$ljoin),
              name=grobname("rect"))
    upViewport(depth)
}
