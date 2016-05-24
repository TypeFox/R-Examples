
setUpUsr <- function(usr) {
    windowvp <- viewport(xscale=usr[1:2], yscale=usr[3:4],
                         name=vpname("window"))
    downViewport(vpPath(vpname("root"), vpname("inner"), vpname("figure"),
                        vpname("plot")), strict=TRUE)
    pushViewport(windowvp)
    upViewport(2)
    downViewport(vpname("plot", clip=TRUE), strict=TRUE)
    pushViewport(windowvp)
    upViewport(3)
    downViewport(vpPath(vpname("figure", clip=TRUE), vpname("plot")),
                 strict=TRUE)
    pushViewport(windowvp)
    upViewport(5)
}

# C_plot_window(xlim, ylim, log, asp, ...)
C_plot_window <- function(x) {
    dev.set(recordDev())
    # NOTE: This takes care of 'asp' (and 'log') by setting par("usr")
    #       appropriately and then we just cream those settings off
    #       for the 'grid' viewports
    do.call("plot.window", x[-1])
    usr <- par("usr")
    dev.set(playDev())
    incrementWindowIndex()
    # Align windowPlotAlpha with plotAlpha
    setWindowPlotAlpha(plotAlpha())
    setUpUsr(usr)
}

