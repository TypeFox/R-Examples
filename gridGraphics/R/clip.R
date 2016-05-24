
# C_clip(x1, x2, y1, y2)

# Just record this clipping setting and enforce it whenever subsequently
# descend into window viewport
C_clip <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:5)])
    dev.set(playDev())
    x1 <- tx(x[[2]], par)
    x2 <- tx(x[[3]], par)
    y1 <- ty(x[[4]], par)
    y2 <- ty(x[[5]], par)
    setClip(x1, y1, x2 - x1, y2 - y1)
}

# Navigate to the correct viewport based on 'xpd' setting
# End up in either "plot" or "window" viewport

gotovp <- function(xpd, end="window") {
    root <- vpname("root")
    inner <- vpname("inner")
    if (is.na(xpd)) {
        figure <- vpname("figure")
        plot <- vpname("plot")
        window <- vpname("window")
        windowplot <- vpname("windowplot")
    } else if (xpd) {
        figure <- vpname("figure", clip=TRUE)
        plot <- vpname("plot")
        window <- vpname("window")
        windowplot <- vpname("windowplot")
    } else {
        figure <- vpname("figure")
        plot <- vpname("plot", clip=TRUE)
        window <- vpname("window")
        windowplot <- vpname("windowplot", clip=TRUE)
    }
    # NOTE that the "window" vp goes via a separate "window" "plot" vp
    #      so that box() can go to one "plot" vp and text() et al can
    #      go to a different "plot" vp (e.g., following a par(mar))
    path <- switch(end,
                   window=vpPath(root, inner, figure, windowplot, window),
                   plot=vpPath(root, inner, figure, plot),
                   figure=vpPath(root, inner, figure),
                   inner=vpPath(root, inner),
                   outer=vpPath(root))
    depth <- downViewport(path, strict=TRUE)
    if (end == "window" && !is.null(clipRegion <- getClip())) {
        grid.clip(clipRegion[1], clipRegion[2], clipRegion[3], clipRegion[4],
                  default.units="native", just=c("left", "bottom"),
                  name=grobname("clip"))
    }
    depth
}


