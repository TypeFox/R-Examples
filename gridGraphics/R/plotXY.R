
# C_plotXY(xy, type, pch, lty, col, bg, cex, lwd, ...)
C_plotXY <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:9)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    xx <- tx(x[[2]]$x, par)
    yy <- ty(x[[2]]$y, par)
    type <- x[[3]]
    pch <- FixupPch(x[[4]], par$pch)
    lty <- FixupLty(x[[5]], par$lty)
    col <- FixupCol(x[[6]], 0, par$bg)
    bg <- FixupCol(x[[7]], NA, par$bg)
    # NOTE: cex multiplied by "base" cex
    cex <- FixupCex(x[[8]]*par$cex, 1)
    lwd <- FixupLwd(x[[9]], par$lwd)
    switch(type,
           n={ }, # do nothing
           p=points(xx, yy, pch, col, bg, cex, lwd, par),
           l=lines(xx, yy, lty, col, lwd, par),
           s=step(xx, yy, lty, col, lwd, par),
           S=Step(xx, yy, lty, col, lwd, par),
           h=bar(xx, yy, lty, col, lwd, par),
           c=brokenlines(xx, yy, lty, col, lwd, par),
           o={ lines(xx, yy, lty, col, lwd, par);
               points(xx, yy, pch, col, bg, cex, lwd, par) },
           b={ brokenlines(xx, yy, lty, col, lwd, par);
               points(xx, yy, pch, col, bg, cex, lwd, par) })
    upViewport(depth)
}

points <- function(x, y, pch, col, bg, cex, lwd, par) {
    grid.points(x, y, default.units="native",
                #  GSTR_0  dpptr(dd)->scale * dd->dev->cra[1] * 0.5 * dd->dev->ipr[1] * gpptr(dd)->cex
                size=unit(par$cin[2]*0.5*cex, "in"), pch=pch,
                gp=gpar(lty="solid", col=col, fill=bg, lwd=lwd, cex=cex,
                    fontface=par$font),
                name=grobname("points"))
}

lines <- function(x, y, lty, col, lwd, par) {
    grid.lines(x, y, default.units="native",
               gp=gpar(lty=lty, col=col, lwd=lwd,
                   lineend=par$lend, linemitre=par$lmitre, linejoin=par$ljoin),
               name=grobname("lines"))
}

step <- function(x, y, lty, col, lwd, par) {
    n <- length(x)
    grid.lines(rep(x, each=2)[-1],
               rep(y, each=2, length.out=2*n - 1),
               default.units="native",
               gp=gpar(lty=lty, col=col, lwd=lwd,
                   lineend=par$lend, linemitre=par$lmitre, linejoin=par$ljoin),
               name=grobname("step"))
}

Step <- function(x, y, lty, col, lwd, par) {
    n <- length(x)
    grid.lines(rep(x, each=2, length.out=2*n - 1),
               rep(y, each=2)[-1],
               default.units="native",
               gp=gpar(lty=lty, col=col, lwd=lwd,
                   lineend=par$lend, linemitre=par$lmitre, linejoin=par$ljoin),
               name=grobname("Step"))
}

bar <- function(x, y, lty, col, lwd, par) {
    if (par$ylog) {
        root <- par$usr[3]
    } else {
        root <- 0
    }
    grid.segments(x, root, x, y, default.units="native",
                  gp=gpar(lty=lty, col=col, lwd=lwd,
                      lineend=par$lend, linemitre=par$lmitre,
                      linejoin=par$ljoin),
                  name=grobname("spike"))
}

brokenlines <- function(x, y, lty, col, lwd, par) {
    d <- 0.5*par$cin[2]*par$cex
    xx <- convertX(unit(x, "native"), "in", valueOnly=TRUE)
    yy <- convertY(unit(y, "native"), "in", valueOnly=TRUE)
    dx <- diff(xx)
    dy <- diff(yy)
    hypot <- sqrt(dx^2 + dy^2)
    # If not enough room, setting to NA will mean no segment drawn
    f <- ifelse(d < 0.5*hypot, d/hypot, NA)
    n <- length(x)
    sx <- xx[-n] + f*dx
    sy <- yy[-n] + f*dy
    ex <- xx[-1] - f*dx
    ey <- yy[-1] - f*dy
    grid.segments(sx, sy, ex, ey, 
                  default.units="in",
                  gp=gpar(lty=lty, col=col, lwd=lwd,
                      lineend=par$lend, linemitre=par$lmitre,
                      linejoin=par$ljoin),
                  name=grobname("brokenline"))
}

