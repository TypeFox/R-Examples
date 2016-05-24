
###################################################
### code chunk number 59: Controls.Rnw:390-396
###################################################
window <- gwindow("gedit example", visible = FALSE) 
group <- ggroup(cont = window)
glabel("State name:", cont = group)
entry <- gedit("", cont = group)
entry[] <- state.name
visible(window) <- TRUE

