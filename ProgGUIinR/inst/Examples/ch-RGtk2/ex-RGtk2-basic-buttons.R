

###################################################
### code chunk number 93: ButtonConstructors
###################################################
window <- gtkWindow(show = FALSE)
window$setTitle("Various buttons")
window$setDefaultSize(400, 25)
hbox <- gtkHBox(homogeneous = FALSE, spacing = 5)
window$add(hbox)
button <- gtkButtonNew() 
button$setLabel("long way")
hbox$packStart(button)
hbox$packStart(gtkButton(label = "label only") )
hbox$packStart(gtkButton(stock.id = "gtk-ok") )
hbox$packStart(gtkButtonNewWithMnemonic("_Mnemonic") )
window$show()
