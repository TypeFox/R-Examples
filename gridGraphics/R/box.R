
# C_box(which, lty, ...)
C_box <- function(x) {
    dev.set(recordDev())
    # NOTE: although 'lty' is passed in, it is not explicitly handled
    par <- currentPar(x[-(1:2)])
    dev.set(playDev())
    which <- x[[2]]
    if (which == 1) { # "plot"
        depth <- gotovp(NA, "plot")
        # NOTE: copy GBox which draws *polygon* (not rect) AND
        #       explicitly sets fill to NA
        xy <- switch(par$bty,
                     "o"=,
                     "O"=list(x=c(0, 1, 1, 0), y=c(0, 0, 1, 1)),
                     "l"=,
                     "L"=list(x=c(0, 0, 1), y=c(1, 0, 0)),
                     "7"=list(x=c(0, 1, 1), y=c(1, 1, 0)),
                     "c"=,
                     "C"=,
                     "["=list(x=c(1, 0, 0, 1), y=c(1, 1, 0, 0)),
                     "]"=list(x=c(0, 1, 1, 0), y=c(1, 1, 0, 0)),
                     "u"=,
                     "U"=list(x=c(0, 0, 1, 1), y=c(1, 0, 0, 1)))
        if (par$bty %in% c("n", "N")) {
            # do nothing
        } else if (par$bty %in% c("o", "O")) {
            grid.polygon(xy$x, xy$y,
                         gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd, fill=NA,
                             lineend=par$lend, linemitre=par$lmitre,
                             linejoin=par$ljoin),
                         name=grobname("box"))
        } else {
            grid.lines(xy$x, xy$y,
                       gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd,
                           lineend=par$lend, linemitre=par$lmitre,
                           linejoin=par$ljoin),
                       name=grobname("box"))
        }
    } else if (which == 2) { # "figure"
        depth <- gotovp(NA, "figure")
        grid.polygon(c(0, 1, 1, 0), c(0, 0, 1, 1),
                     gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd, fill=NA,
                         lineend=par$lend, linemitre=par$lmitre,
                         linejoin=par$ljoin),
                     name=grobname("box-figure"))
    } else if (which == 3) { # "inner"
        depth <- gotovp(NA, "inner")
        grid.polygon(c(0, 1, 1, 0), c(0, 0, 1, 1),
                     gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd, fill=NA,
                         lineend=par$lend, linemitre=par$lmitre,
                         linejoin=par$ljoin),
                     name=grobname("box-inner"))
    } else { # "outer"
        depth <- gotovp(NA, "outer")
        grid.polygon(c(0, 1, 1, 0), c(0, 0, 1, 1),
                     gp=gpar(col=par$col, lty=par$lty, lwd=par$lwd, fill=NA,
                         lineend=par$lend, linemitre=par$lmitre,
                         linejoin=par$ljoin),
                     name=grobname("box-outer"))        
    }
    upViewport(depth)
}

