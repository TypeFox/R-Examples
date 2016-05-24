### R code from vignette source 'ex-RGtk2-basic-menu.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-basic-menu.Rnw:5-8
###################################################
mb <- gtkMenuBar()
fileMi <- gtkMenuItemNewWithMnemonic(label="_File")
mb$append(fileMi)


###################################################
### code chunk number 2: ex-RGtk2-basic-menu.Rnw:12-14
###################################################
fileMb <- gtkMenu()
fileMi$setSubmenu(fileMb)


###################################################
### code chunk number 3: ex-RGtk2-basic-menu.Rnw:17-18
###################################################
open <- gtkMenuItem(label="open")


###################################################
### code chunk number 4: ex-RGtk2-basic-menu.Rnw:22-24
###################################################
saveAction <- gtkAction("save", "save", "Save object", "gtk-save")
save <- saveAction$CreateMenuItem()


###################################################
### code chunk number 5: ex-RGtk2-basic-menu.Rnw:28-31
###################################################
quit <- gtkImageMenuItem(label="quit")
quit$setImage(gtkImageNewFromStock("gtk-quit", 
              size=GtkIconSize["menu"]))


###################################################
### code chunk number 6: ex-RGtk2-basic-menu.Rnw:35-37
###################################################
happy <- gtkCheckMenuItem(label="happy")
happy$setActive(TRUE)


###################################################
### code chunk number 7: ex-RGtk2-basic-menu.Rnw:41-45
###################################################
items <- list(open, save, happy, gtkSeparatorMenuItem(), quit)
Qt <- sapply(items, function(i) {
       fileMb$append(i)
     })


###################################################
### code chunk number 8: ex-RGtk2-basic-menu.Rnw:48-54
###################################################
ID <- gSignalConnect(happy, "toggled", function(b,data) {
  if(b$getActive())
    print("User is now happy ;)")
  else
    print("User is unhappy ;(")
})


###################################################
### code chunk number 9: ex-RGtk2-basic-menu.Rnw:57-62
###################################################
QT <- sapply(list(open, quit, saveAction), function(i) 
       gSignalConnect(i, "activate", f=function(mi, data) {
         cat("item selected\n")
       })
       )


###################################################
### code chunk number 10: makeGUI
###################################################
## We make as simple GUI for the menubar.
 w <- gtkWindow(show=FALSE)
w['title'] <- "Menubar example"
w$add(mb)
w$ShowAll()


