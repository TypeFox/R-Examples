
###################################################
### code chunk number 61: BasicContainers.Rnw:119-125
###################################################
## set up window etc. for example
library(tcltk)
window <- tktoplevel()
tkwm.title(window, "Simple text entry example")
entry <- tktext(window)
tkpack(entry)


###################################################
### code chunk number 62: example-wm-delete-window
###################################################
tkwm.protocol(window, "WM_DELETE_WINDOW", function() {
  modified <- tcl(entry, "edit", "modified")
  if(as.logical(modified)) {
    response <- 
      tkmessageBox(icon = "question",
                   message = "Really close?",
                   detail = "Changes need to be saved",
                   type = "yesno",
                   parent = window)
    if(as.character(response) == "no")
      return()
  }
  tkdestroy(window)                     # otherwise close
})
