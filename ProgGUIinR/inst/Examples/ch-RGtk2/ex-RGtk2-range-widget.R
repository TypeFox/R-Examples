### R code from vignette source 'ex-RGtk2-range-widget.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-range-widget.Rnw:15-17
###################################################
## make a range widget combining both a slider and spinbutton to choose a number
library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-range-widget.Rnw:22-23
###################################################
from <- 0; to <- 100; by <- 1


###################################################
### code chunk number 3: ex-RGtk2-range-widget.Rnw:30-34
###################################################
slider <- gtkHScale(min = from, max = to, step = by)
slider['draw-value'] <- FALSE
adjustment <- slider$getAdjustment()
spinbutton <- gtkSpinButton(adjustment = adjustment)


###################################################
### code chunk number 4: ex-RGtk2-range-widget.Rnw:41-44
###################################################
hbox <- gtkHBox()
hbox$packStart(slider, expand = TRUE, fill = TRUE, padding = 5)
hbox$packStart(spinbutton, expand = FALSE, padding = 5)


###################################################
### code chunk number 5: ex-RGtk2-range-widget.Rnw:48-53
###################################################
w <- gtkWindow(show=FALSE)
w['title'] <- "Example of a range widget"
w$setSizeRequest(width=200, height=-1)
w$add(hbox)
w$show()


