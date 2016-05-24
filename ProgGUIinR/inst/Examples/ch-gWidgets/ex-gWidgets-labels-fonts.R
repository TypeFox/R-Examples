###################################################
### code chunk number 52: Controls.Rnw:218-228
###################################################
window <- gwindow("label example")
frame <- gframe("Summary statistics:", cont = window)
lyt <- glayout(cont = frame)
lyt[1,1] <- glabel("xbar:", cont = lyt)
lyt[1,2] <- gedit("", cont = lyt)
lyt[2,1] <- glabel("s:", cont = lyt)
lyt[2,2] <- gedit("", cont = lyt)
sapply(lyt[,1], function(i) {
  font(i) <- c(weight = "bold", color = "blue")
})
