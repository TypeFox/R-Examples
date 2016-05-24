
C_plot_new <- function(x) {
    # recordDev
    dev.set(recordDev())
    page <- par("page")
    plot.new()
    par <- par()
    # playDev
    dev.set(playDev())
    initFigureAlpha()
    incrementPlotIndex()
    initPlotAlpha()
    initWindowIndex()
    initWindowAlpha()
    initWindowPlotAlpha()
    initClip()
    nvp <- 0
    if (page) {
        if (get("newpage", .gridGraphicsEnv))
            grid.newpage()
        # If there is a non-transparent par(bg) in effect, we need
        # to draw an opaque background rect
        if (par$bg != "transparent") {
            grid.rect(width=1.5, height=1.5,
                      gp=gpar(col=NA, fill=par$bg),
                      name=grobname("background", unique=TRUE))
        }
        pushViewport(viewport(# gp=gparFromPar(par[gparParNames]),
                              name=vpname("root")))
        upViewport()
        setUpInner(par)
    } else {
        setUpFigure(par)
    }
}

setUpInner <- function(par) {
    omi <- par$omi
    innervp <- viewport(x=unit(omi[2], "inches"),
                        y=unit(omi[1], "inches"),
                        width=unit(1, "npc") -
                        unit(omi[2], "inches") -
                        unit(omi[4], "inches"),
                        height=unit(1, "npc") -
                        unit(omi[1], "inches") -
                        unit(omi[3], "inches"),
                        just=c("left", "bottom"),
                        name=vpname("inner"))
    downViewport(vpname("root"), strict=TRUE)
    pushViewport(innervp)
    upViewport(2)
    setUpFigure(par)
}

setUpFigure <- function(par) {
    fig <- par$fig
    figurevp <- viewport(x=unit(fig[1], "npc"),
                         y=unit(fig[3], "npc"),
                         width=unit(fig[2] - fig[1], "npc"),
                         height=unit(fig[4] - fig[3], "npc"),
                         just=c("left", "bottom"),
                         name=vpname("figure"))
    figurevpclip <- viewport(x=unit(fig[1], "npc"),
                             y=unit(fig[3], "npc"),
                             width=unit(fig[2] - fig[1], "npc"),
                             height=unit(fig[4] - fig[3], "npc"),
                             just=c("left", "bottom"),
                             clip=TRUE,
                             name=vpname("figure", clip=TRUE))    
    downViewport(vpPath(vpname("root"), vpname("inner")), strict=TRUE)
    pushViewport(figurevp)
    upViewport()
    pushViewport(figurevpclip)
    upViewport(3)
    setUpPlot(par)
}

setUpPlot <- function(par) {
    plt <- par$plt
    plotvp <- viewport(x=unit(plt[1], "npc"),
                       y=unit(plt[3], "npc"),
                       width=unit(plt[2] - plt[1], "npc"),
                       height=unit(plt[4] - plt[3], "npc"),
                       just=c("left", "bottom"),
                       name=vpname("plot"))
    plotvpclip <- viewport(x=unit(plt[1], "npc"),
                           y=unit(plt[3], "npc"),
                           width=unit(plt[2] - plt[1], "npc"),
                           height=unit(plt[4] - plt[3], "npc"),
                           just=c("left", "bottom"),
                           clip=TRUE,
                           name=vpname("plot", clip=TRUE))
    downViewport(vpPath(vpname("root"), vpname("inner"), vpname("figure")),
                 strict=TRUE)
    pushViewport(plotvp)
    upViewport(1)
    pushViewport(plotvpclip)
    upViewport(2)
    downViewport(vpname("figure", clip=TRUE), strict=TRUE)
    pushViewport(plotvp)
    upViewport(4)
    setUpUsr(par$usr)
}
