
# C_abline(a, b, h, v, untf, col, lty, lwd, ...)
C_abline <- function(x) {
    # TODO: handle 'untf'
    dev.set(recordDev())
    par <- currentPar(x[-(1:9)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    a <- x[[2]]
    b <- x[[3]]
    h <- ty(x[[4]], par)
    v <- tx(x[[5]], par)
    untf <- x[[6]]
    col <- FixupCol(x[[7]], NA, par$bg)
    lty <- FixupLty(x[[8]], par$lty)
    lwd <- FixupLwd(x[[9]], par$lwd)
    if (!is.null(a)) {
        if (is.null(b)) {
            a <- a[1]
            b <- a[2]
        }
        xx <- par$usr[1:2]
        if (untf && (par$xlog || par$ylog)) {
            # Draw a curve instead of straight line
            NS <- 100
            if (par$xlog) {
                xx <- 10^xx
                xf <- xx[2]/.Machine$double.xmax
                xxx <- xx[1] <- max(xx[1], 1.01*xf)
                xf <- (xx[2]/xx[1])^(1/NS)
                xxx <- xxx*xf^(0:(NS - 1))
                xxx[NS + 1] <- xx[2]
            } else {
                xstep <- (xx[2] - xx[1])/NS
                xxx <- seq(xx[1], xx[2], xstep)
            }
            yyy <- b*xxx + a
            if (par$xlog) {
                xxx <- log10(xxx)
                # xxx[xxx <= 0] <- NA
            }
            if (par$ylog) {
                yyy <- log10(yyy)
                # yyy[yyy <= 0] <- NA
            }
            xx <- xxx
            yy <- yyy
        } else {
            # If 'xlog' then par$usr is already "logged"
            # ditto for 'ylog'
            yy <- a + b*xx
        }
        # TODO: this will have to be smarter to handle drawing outside
        #       the plot region
        grid.lines(xx, yy, default.units="native",
                   gp=gpar(col=col, lty=lty, lwd=lwd,
                       lineend=par$lend, linemitre=par$lmitre,
                       linejoin=par$ljoin),
                   name=grobname("abline-ab"))
    }
    if (!is.null(h)) {
        grid.segments(0, unit(h, "native"), 1, unit(h, "native"),
                      gp=gpar(col=col, lty=lty, lwd=lwd,
                          lineend=par$lend, linemitre=par$lmitre,
                          linejoin=par$ljoin),
                      name=grobname("abline-h"))
    }
    if (!is.null(v)) {
        grid.segments(unit(v, "native"), 0, unit(v, "native"), 1, 
                      gp=gpar(col=col, lty=lty, lwd=lwd,
                          lineend=par$lend, linemitre=par$lmitre,
                          linejoin=par$ljoin),
                      name=grobname("abline-v"))
    }
    upViewport(depth)
}
