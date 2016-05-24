### R code from vignette source 'ex-RGtk2-editable-label.Rnw'

###################################################
### code chunk number 1: EditableLabel
###################################################
w <- gtkWindow(); w$setTitle("Editable label")
evb <- gtkEventBox(); 
w$add(evb); 
e <- gtkEntry()
l <- gtkLabel("Click me to edit")
evb$setData("entry", e); evb$setData("label", l)
evb$add(l)

ID <- gSignalConnect(evb, "button-press-event", function(w, e, ...) {
  label <- w$getData("label"); entry <- w$getData("entry")
  entry$setText(label$getText())
  w$remove(label);  w$add(entry) # swap
})
ID <- gSignalConnect(e,"activate", function(userData, w, ...) {
  evb = userData$evb; label <- evb$getData("label")
  label$setText(w$getText())
  evb$remove(w); evb$add(label) #swap
},
                     data=list(evb=evb),
                     user.data.first=TRUE)


