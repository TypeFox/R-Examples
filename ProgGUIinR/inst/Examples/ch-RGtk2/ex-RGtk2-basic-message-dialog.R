###################################################
### code chunk number 61: Pre-defined-dialogs.Rnw:21-29
###################################################
window <- gtkWindow(); window['title'] <- "Parent window"
#
dialog <- gtkMessageDialog(parent=window, 
                           flags="destroy-with-parent",
                           type="question", 
                           buttons="ok",
                           "My message")
dialog['secondary-text'] <- "A secondary message"


###################################################
### code chunk number 62: Pre-defined-dialogs.Rnw:49-58
###################################################
response <- dialog$run()
if(response == GtkResponseType["cancel"] ||
   response == GtkResponseType["close"] ||
   response == GtkResponseType["delete-event"]) {
  ## pass
} else if(response == GtkResponseType["ok"]) {
  message("Ok")
}
dialog$destroy()

