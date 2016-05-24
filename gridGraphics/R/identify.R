
C_identify <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:9)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    index <- x[[2]]
    pos <- x[[3]]
    xx <- unit(x[[4]], "native")
    yy <- unit(x[[5]], "native")
    offset <- unit(x[[6]]*par$cin[2]*par$cex, "in")
    label <- x[[7]]
    draw <- x[[8]]
    for (i in seq_along(xx)) {
        xxx <- switch(pos[i] + 1,
                      xx[i],
                      xx[i], xx[i] - offset, xx[i], xx[i] + offset)
        yyy <- switch(pos[i] + 1,
                      yy[i],
                      yy[i] - offset, yy[i], yy[i] + offset, yy[i])
        xadj <- switch(pos[i] + 1, 0, 0.5, 1, 0.5, 0)
        # 0.3333 comes from dev->yCharOffset
        yadj <- switch(pos[i] + 1, 0, 1 - (0.5 - 0.3333), 0.3333, 0, 0.3333)
        if (index[i] && draw) {
            grid.text(label, xxx, yyy, default.units="native",
                      hjust=xadj, vjust=yadj,
                      gp=gpar(cex=par$cex, col=par$col, fontface=par$font,
                          fontfamily=par$family, lineheight=par$lheight))
        }
    }
    upViewport(depth)
}

