
# C_symbols(x, y, type, data, inches, bg, fg, ...)

C_symbols <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:8)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    xx <- tx(x[[2]], par)
    yy <- ty(x[[3]], par)
    type <- x[[4]]
    p <- x[[5]]
    # FALSE becomes 0, TRUE becomes 1
    inches <- as.numeric(x[[6]])
    if (!is.finite(inches) || inches < 0) {
        inches <- 0
    }
    bg <- FixupCol(x[[7]], NA, par$bg)
    fg <- FixupCol(x[[8]], NA, par$bg)
    if (type == 1) { # circles
        prange <- range(p, na.rm=TRUE)
        if (inches > 0) {
            r <- unit(p*inches/prange[2], "in")
        } else {
            r <- convertWidth(unit(p, "native"), "in")
        }
        grid.circle(xx, yy, r, default.units="native",
                    gp=gpar(col=fg, fill=bg, lty=par$lty, lwd=par$lwd),
                    name=grobname("symbols-circle"))
    } else if (type == 2) { # squares
        prange <- range(p, na.rm=TRUE)
        if (inches > 0) {
            w <- unit(p*inches/prange[2], "in")
        } else {
            w <- convertWidth(unit(p, "native"), "in")
        }
        grid.rect(xx, yy, width=w, height=w, default.units="native",
                  gp=gpar(col=fg, fill=bg, lty=par$lty, lwd=par$lwd,
                      lineend=par$lend, linemitre=par$lmitre,
                      linejoin=par$ljoin),
                  name=grobname("symbols-square"))
    } else if (type == 3) { # rectangles
        prange <- range(p, na.rm=TRUE)
        if (inches > 0) {
            w <- unit(p[, 1]*inches/prange[2], "in")
            h <- unit(p[, 2]*inches/prange[2], "in")
        } else {
            w <- unit(p[, 1], "native")
            h <- unit(p[, 2], "native")
        }
        grid.rect(xx, yy, width=w, height=h, default.units="native",
                  gp=gpar(col=fg, fill=bg, lty=par$lty, lwd=par$lwd,
                      lineend=par$lend, linemitre=par$lmitre,
                      linejoin=par$ljoin),
                  name=grobname("symbols-rect")) 
    } else if (type == 4) { # stars
        prange <- range(p, na.rm=TRUE)
        nc <- ncol(p)
        nr <- nrow(p)
        p[is.na(p)] <- 0
        if (inches > 0) {
            r <- p*inches/prange[2]
        } else {
            r <- convertWidth(unit(p, "native"), "in", valueOnly=TRUE)
        }
        xc <- rep(convertX(unit(xx, "native"), "in", valueOnly=TRUE), nc)
        yc <- rep(convertX(unit(yy, "native"), "in", valueOnly=TRUE), nc)
        t <- seq(0, 2*pi, length.out=nc + 1)[-1]
        grid.polygon(xc + r*cos(t),  yc + r*sin(t), default.units="in",
                     id.lengths=rep(1:nc, each=nr),
                     gp=gpar(col=fg, fill=bg, lty=par$lty, lwd=par$lwd,
                         lineend=par$lend, linemitre=par$lmitre,
                         linejoin=par$ljoin),
                     name=grobname("symbols-star")) 
    } else if (type == 5) { # thermometers
        prange <- range(p[, 1:2], na.rm=TRUE)
        nc <- ncol(p)
        if (nc < 4) {
            p <- cbind(p, 0)
        }
        if (inches > 0) {
            w <- unit(p[, 1]*inches/prange[2], "in")
            h <- unit(p[, 2]*inches/prange[2], "in")
        } else {
            w <- convertWidth(unit(p[, 1], "native"), "in")
            h <- convertHeight(unit(p[, 2], "native"), "in")
        }
        grid.rect(xx, yy, width=w, height=h, default.units="native",
                  gp=gpar(col=fg, fill=bg, lty=par$lty, lwd=par$lwd,
                      lineend=par$lend, linemitre=par$lmitre,
                      linejoin=par$ljoin),
                  name=grobname("symbols-thermo-box")) 
        grid.rect(unit(xx, "native"),
                  unit(yy, "native") - (1 - 2*p[, 3])*0.5*h,
                  width=w,
                  height=(p[, 3] - p[, 4])*h,
                  just="top",
                  gp=gpar(col=fg, fill=fg, lty=par$lty, lwd=par$lwd,
                      lineend=par$lend, linemitre=par$lmitre,
                      linejoin=par$ljoin),
                  name=grobname("symbols-thermo-fill"))
        grid.segments(unit(xx, "native") + 0.5*w,
                      unit(yy, "native"),
                      unit(xx, "native") + 0.75*w,
                      unit(yy, "native"),
                      gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd,
                          lineend=par$lend, linemitre=par$lmitre,
                          linejoin=par$ljoin),
                      name=grobname("symbols-thermo-whisker-right"))
        grid.segments(unit(xx, "native") - 0.5*w,
                      unit(yy, "native"),
                      unit(xx, "native") - 0.75*w,
                      unit(yy, "native"),
                      gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd,
                          lineend=par$lend, linemitre=par$lmitre,
                          linejoin=par$ljoin),
                      name=grobname("symbols-thermo-whisker-left"))       
    } else if (type == 6) { # boxplots
        prange <- range(p[, 1:4], na.rm=TRUE)
        if (inches > 0) {
            w <- unit(p[, 1]*inches/prange[2], "in")
            h <- unit(p[, 2]*inches/prange[2], "in")
            lw <- unit(p[, 3]*inches/prange[2], "in")
            uw <- unit(p[, 4]*inches/prange[2], "in")
        } else {
            w <- unit(p[, 1], "native")
            h <- unit(p[, 2], "native")
            lw <- convertHeight(unit(p[, 3], "native"), "in")
            uw <- convertHeight(unit(p[, 4], "native"), "in")
        }
        grid.rect(xx, yy, width=w, height=h, default.units="native",
                  gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd,
                      lineend=par$lend, linemitre=par$lmitre,
                      linejoin=par$ljoin),
                  name=grobname("symbols-boxplot-box"))       
        grid.segments(unit(xx, "native"),
                      unit(yy, "native") - 0.5*h,
                      unit(xx, "native"),
                      unit(yy, "native") - 0.5*h - lw,
                      gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd,
                          lineend=par$lend, linemitre=par$lmitre,
                          linejoin=par$ljoin),
                      name=grobname("symbols-boxplot-lower-whisker"))       
        grid.segments(unit(xx, "native"),
                      unit(yy, "native") + 0.5*h,
                      unit(xx, "native"),
                      unit(yy, "native") + 0.5*h + uw,
                      gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd,
                          lineend=par$lend, linemitre=par$lmitre,
                          linejoin=par$ljoin),
                      name=grobname("symbols-boxplot-upper-whisker"))       
        grid.segments(unit(xx, "native") - 0.5*w,
                      unit(yy, "native") - (1 - 2*p[, 3])*0.5*h,
                      unit(xx, "native") + 0.5*w,
                      unit(yy, "native") - (1 - 2*p[, 3])*0.5*h,
                      gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd,
                          lineend=par$lend, linemitre=par$lmitre,
                          linejoin=par$ljoin),
                      name=grobname("symbols-boxplot-median"))
    }
    upViewport(depth)
}
