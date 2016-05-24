###################################################
### code chunk number 84: slider
###################################################
window <- gwindow("Add a label to the slider", visible=FALSE)
group <- ggroup(cont = window, expand = TRUE)
slider <- gslider(from = 0, to = 100, by = 1, cont = group, 
                  expand = TRUE)
label <- glabel(sprintf("%3d", svalue(slider)), cont = group)
font(label) <- c(family = "monospace")
addHandlerChanged(slider, function(h,...) {
  svalue(h$action) <- sprintf("%3d", svalue(h$obj))
  }, action = label)
visible(window) <- TRUE
