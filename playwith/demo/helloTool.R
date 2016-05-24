
library(playwith)

## 1. A simple action to show a "Hello world" dialog box.

hello_handler <- function(widget, playState)
    gmessage("Hello world.")

helloTool <- list("Hello", label = "Say Hello",
                  callback = hello_handler)

playwith(plot(1:10), tools = list(helloTool))


## 2. A toggle button to draw "G'day universe" on the plot.

## This "callback" is run when the action is activated
## (from toolbar, menu, or keyboard shortcut).
gday_handler <- function(widget, playState) {
    if (widget["active"]) {
        ## tool was turned on: draw the text
        playState$tmp$gdayActive <- TRUE
        drawGday(playState)
    } else {
        ## turned off; re-draw plot to remove text
        playState$tmp$gdayActive <- FALSE
        playReplot(playState)
    }
}

## This is an "update action", called after plotting.
drawGday <- function(playState) {
    if (isTRUE(playState$tmp$gdayActive)) {
        ## draw text centred on the page
        grid.text("G'day universe.", gp = gpar(cex=2))
    }
}

gdayTool <- list("Gday", "gtk-yes", "Draw G'day", "F5",
                 "Overlay text on the plot", gday_handler, FALSE,
                 update.action = drawGday)

playwith(plot(1:10), tools = list(gdayTool))
## Note the toolbar button (see playwith.options("custom.toolbar"))
## and item in the Tools menu (with keyboard shortcut F5).
