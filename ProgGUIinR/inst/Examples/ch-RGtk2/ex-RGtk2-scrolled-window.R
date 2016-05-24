### R code from vignette source 'ex-RGtk2-scrolled-window.Rnw'

###################################################
### code chunk number 1: ScrolledWindowExample
###################################################
## Simplistic example of a scrolled window
## there are better uses.
library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-scrolled-window.Rnw:17-23
###################################################
g <- gtkVBox(spacing=0)
sapply(state.name, function(i) {
  l <- gtkLabel(i)
  l['xalign'] <- 0; l['xpad'] <- 10
  g$packStart(l, expand=TRUE, fill=TRUE)
})


###################################################
### code chunk number 3: ex-RGtk2-scrolled-window.Rnw:28-31
###################################################
sw <- gtkScrolledWindow()
sw$setPolicy("never","automatic")
sw$addWithViewport(g)          # just "Add" for text, tree, ...


###################################################
### code chunk number 4: ex-RGtk2-scrolled-window.Rnw:34-39
###################################################
w <- gtkWindow(show=FALSE)
w$setTitle("Scrolled window example")
w$setSizeRequest(-1, 300)
w$add(sw)
w$show()


