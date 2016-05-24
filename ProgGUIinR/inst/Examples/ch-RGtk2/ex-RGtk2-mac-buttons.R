### R code from vignette source 'ex-RGtk2-mac-buttons.Rnw'

###################################################
### code chunk number 1: MacOSXstyleButton
###################################################
## not shown
window <- gtkWindow(show=FALSE)
window$setTitle("MAC OS X style buttons")
fg <- gtkVBox()
fg$setSizeRequest(width=800, height=-1)
window$add(fg)

hbox <- gtkHBox()
fg$packStart(hbox, padding=15)              # for size grip


###################################################
### code chunk number 2: StockButtons
###################################################
ok <- gtkButton(stock.id="gtk-ok")
cancel <- gtkButton(stock.id="gtk-cancel")
delete <- gtkButton(stock.id="gtk-delete")


###################################################
### code chunk number 3: macButtonPack
###################################################
hbox$packEnd(ok, padding = 0)
hbox$packEnd(cancel, padding = 12)
hbox$packEnd(delete, padding = 12)
hbox$packEnd(gtkLabel(""), expand = TRUE, fill = TRUE)


###################################################
### code chunk number 4: ex-RGtk2-mac-buttons.Rnw:56-57
###################################################
ok$grabFocus()


###################################################
### code chunk number 5: ex-RGtk2-mac-buttons.Rnw:60-62
###################################################
## not shown
window$showAll()


###################################################
### code chunk number 6: gtkHButtonBoxExample
###################################################
## not shown
## Had we only wanted to use a button box
button_box <- gtkHButtonBox()
fg$packStart(button_box, padding=15)              # for size grip

button_box$add(gtkButton(stock.id="gtk-delete"))
button_box$add(gtkButton(stock.id="gtk-cancel"))
button_box$add(gtkButton(stock.id="gtk-ok"))


