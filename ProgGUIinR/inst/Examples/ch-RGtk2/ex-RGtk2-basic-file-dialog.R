
###################################################
### code chunk number 66: openFileDialog
###################################################
dialog <- gtkFileChooserDialog(title = "Open a file", 
                     parent = NULL, action = "open",
                     "gtk-ok", GtkResponseType["ok"],
                     "gtk-cancel", GtkResponseType["cancel"],
                     show = FALSE)


###################################################
### code chunk number 67: Pre-defined-dialogs.Rnw:165-173
###################################################
gSignalConnect(dialog, "response", 
               f = function(dialog, response, data) {
                 if(response == GtkResponseType["ok"]) {
                   filename <- dialog$getFilename()
                   print(filename)
                 }
                 dialog$destroy()
               })


###################################################
### code chunk number 68: Pre-defined-dialogs.Rnw:183-188
###################################################
fileFilter <- gtkFileFilter()
fileFilter$setName("R files")
fileFilter$addPattern("*.R")
fileFilter$addPattern("*.Rdata")
dialog$addFilter(fileFilter)

dialog$run()
