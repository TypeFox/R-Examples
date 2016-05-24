
###################################################
### code chunk number 110: panedWindowExample
###################################################
window <- tktoplevel()
tkwm.title(window, "Paned window example")
paned <- ttkpanedwindow(window, orient = "horizontal")
tkpack(paned, expand = TRUE, fill = "both")
left <- ttklabel(paned, text = "left")
right <- ttklabel(paned, text = "right")
#
tkadd(paned, left, weight = 1)
tkadd(paned, right, weight = 2)


###################################################
### code chunk number 111: BasicContainers.Rnw:732-734
###################################################
tcl(paned, "forget", right)
tkadd(paned, right, weight = 2) ## how to add back


###################################################
### code chunk number 112: BasicContainers.Rnw:743-745
###################################################
width <- as.integer(tkwinfo("width", paned))  # or "height"
tcl(paned, "sashpos", 0, floor(0.75*width))
