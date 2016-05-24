
library(playwith)

## A button to interactively add or remove data points.

points_handler <- function(widget, playState) {
    repeat {
        foo <- playSelectData(playState,
                              prompt=paste(
                              "Click to add a point.",
                              "Ctrl-click to delete.",
                              "Right-click to stop."))
        if (is.null(foo)) return()
        xy <- xyData(playState)
        if (foo$modifiers & RGtk2::GdkModifierType["control-mask"]) {
            ## Ctrl-click: delete data points
            xy$x[foo$which] <- NA
            xy$y[foo$which] <- NA
        } else {
            ## add data point at click location
            xy$x <- c(xy$x, foo$coords$x[1])
            xy$y <- c(xy$y, foo$coords$y[1])
        }
        ## store in local environment
        playState$env$localxy <- xy
        if (playState$is.lattice) {
            ## lattice plot: use formula
            callArg(playState, 1) <- quote(y ~ x)
            callArg(playState, "data") <- quote(localxy)
        } else {
            ## otherwise set first argument to plot
            callArg(playState, 1) <- quote(localxy)
            callArg(playState, "y") <- NULL
        }
        playReplot(playState)
    }
}

pointsTool <- list("Points", "gtk-add", "Add points",
                   callback = points_handler)

ydata <- c(1:4, 2:1, 5:8)
playwith(xyplot(ydata ~ 1:10, type=c("p", "smooth"), pch=8),
	tools = list(pointsTool))
