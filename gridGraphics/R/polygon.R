
# C_polygon(x, y, col, border, lty, ...)

C_polygon <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:6)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    xx <- tx(x[[2]], par)
    yy <- ty(x[[3]], par)
    col <- FixupCol(x[[4]], NA, par$bg)
    border <- FixupCol(x[[5]], par$fg, par$bg)
    lty <- FixupLty(x[[6]], par$lty)
    lty <- ifelse(is.na(lty), par$lty, lty)
    # NOTE: allow for NA values in x/y
    breaks <- which(is.na(xx) | is.na(yy))
    if (length(breaks) == 0) { # Only one polygon
        grid.polygon(xx, yy, default.units="native",
                     gp=gpar(col=border, fill=col, lty=lty, lwd=par$lwd,
                         lineend=par$lend, linemitre=par$lmitre,
                         linejoin=par$ljoin),
                     name=grobname("polygon"))
    } else {
        nb <- length(breaks)
        lengths <- c(breaks[1] - 1,
                     diff(breaks) - 1,
                     length(xx) - breaks[nb])
        grid.polygon(xx[-breaks], yy[-breaks], default.units="native",
                     id.lengths=lengths,
                     gp=gpar(col=border, fill=col, lty=lty, lwd=par$lwd,
                         lineend=par$lend, linemitre=par$lmitre,
                         linejoin=par$ljoin),
                     name=grobname("polygon"))
    }
    upViewport(depth)
}

