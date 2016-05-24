

###################################################
### code chunk number 171: tkspinbox
###################################################
tkspinbox <- function(parent, ...) 
    tkwidget(parent, "tk::spinbox", ...)


###################################################
### code chunk number 172: Widgets.Rnw:552-554
###################################################
window <- tktoplevel()
spin_box <- tkspinbox(window, values = state.name, wrap=TRUE)


###################################################
### code chunk number 173: Widgets.Rnw:558-559
###################################################
spin_box1 <- tkspinbox(window, from=1, to = 10, increment = 1)


###################################################
### code chunk number 174: Widgets.Rnw:562-564
###################################################
tkpack(spin_box)
tkpack(spin_box1)
