
# C_path(x, y, lengths, rule, col, border, lty, ...)

C_path <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:8)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    xx <- tx(x[[2]], par)
    yy <- ty(x[[3]], par)
    lengths <- x[[4]]
    rule <- x[[5]]
    col <- FixupCol(x[[6]], NA, par$bg)
    border <- FixupCol(x[[7]], par$fg, par$bg)
    lty <- FixupLty(x[[8]], par$lty)
    lty <- ifelse(is.na(lty), par$lty, lty)
    grid.path(xx, yy, default="native",
              id.lengths=lengths, rule=rule,
              gp=gpar(col=border, fill=col, lty=lty, lwd=par$lwd,
                  lineend=par$lend, linemitre=par$lmitre,
                  linejoin=par$ljoin),
              name=grobname("path"))
    upViewport(depth)    
}
