
###################################################
### code chunk number 123: ColorSelection
###################################################
window <- tktoplevel()
tkwm.title(window, "Select a color")
frame <- ttkframe(window, padding = c(3,3,3,12))
tkpack(frame, expand = TRUE, fill = "both")
color_well <- tkcanvas(frame, width = 40, height = 16, 
                      background = "#ee11aa",
                      highlightbackground = "#ababab") 
tkpack(color_well)
tkpack(ttklabel(frame, text = "Click color to change"))
#
tkbind(color_well,"<Button-1>", function(W) {
  color <- tcl("tk_chooseColor", parent = W, 
               title = "Set box color")
  color <- tclvalue(color)
  print(color)
  if(nchar(color))
    tkconfigure(W, background = color)
})
