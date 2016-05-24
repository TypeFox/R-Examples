
###################################################
### code chunk number 63: Pre-defined-dialogs.Rnw:86-91
###################################################
dialog <- gtkDialogNewWithButtons(title = "Enter a value", 
                   parent = NULL, flags = 0,
                   "gtk-ok", GtkResponseType["ok"],
                   "gtk-cancel", GtkResponseType["cancel"],
                   show = FALSE)


###################################################
### code chunk number 64: OurDialogsLayout
###################################################
hbox <- gtkHBox()
hbox['spacing'] <- 10
#
hbox$packStart(gtkLabel("Enter a value:"))
entry <- gtkEntry()
hbox$packStart(entry)
#
vbox <- dialog$getContentArea()
vbox$packStart(hbox)


###################################################
### code chunk number 65: connectResponse
###################################################
gSignalConnect(dialog, "response", 
               f=function(dialog, response, user.data) {
                 if(response == GtkResponseType["ok"])
                   print(entry$getText()) # Replace this
                 dialog$Destroy()
               })
dialog$showAll()
dialog$setModal(TRUE)
