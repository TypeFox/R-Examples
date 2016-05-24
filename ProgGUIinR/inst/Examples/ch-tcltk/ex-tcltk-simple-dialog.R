### R code from vignette source 'ex-tcltk-simple-dialog.Rnw'

###################################################
### code chunk number 1: ex-tcltk-simple-dialog.Rnw:11-15
###################################################
## example of  simple non-modal dialog
## much taken from msgbox.tcl in tk source

library(tcltk)


###################################################
### code chunk number 2: ex-tcltk-simple-dialog.Rnw:19-25
###################################################
title <- "message dialog"
message <- "Do you like tcltk so far?"
parent <- NULL
tkimage.create("photo", "::img::tclLogo", 
               file = system.file("images","tclp.gif",
                 package = "ProgGUIinR"))


###################################################
### code chunk number 3: ex-tcltk-simple-dialog.Rnw:30-35
###################################################
window <- tktoplevel()
tkwm.title(window, title)
tkwm.state(window, "withdrawn")
frame <- ttkframe(window,  padding = c(3, 3, 12, 12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 4: ex-tcltk-simple-dialog.Rnw:45-58
###################################################
if(!is.null(parent)) {
  parent_window <- tkwinfo("toplevel", parent)
  if(as.logical(tkwinfo("viewable", parent_window))) {
    tkwm.transient(window, parent)
    ## have fun with OS X
    if(as.character(tcl("tk", "windowingsystem")) == "aqua") {
      tcl("wm","attributes", parent_window, notify = TRUE)
      tkbind(parent_window, "<Enter>", function() 
             tcl("wm","attributes", parent_window, 
                 notify = FALSE))       # stop bounce
    }
  }
}


###################################################
### code chunk number 5: ex-tcltk-simple-dialog.Rnw:63-65
###################################################
label <- ttklabel(frame, image = "::img::tclLogo", padding=5) 
tkpack(label, side = "left")


###################################################
### code chunk number 6: ex-tcltk-simple-dialog.Rnw:71-76
###################################################
frame_1 <- ttkframe(frame)
tkpack(frame_1, expand = TRUE, fill = "both")
#
m <- ttklabel(frame_1, text = message)
tkpack(m, expand = TRUE, fill = "both")


###################################################
### code chunk number 7: ex-tcltk-simple-dialog.Rnw:80-82
###################################################
frame_2 <- ttkframe(frame_1)
tkpack(frame_2)


###################################################
### code chunk number 8: ex-tcltk-simple-dialog.Rnw:86-97
###################################################
ok_callback <- function() {
  print("That's great")
  tkdestroy(w)
}
ok_button <- ttkbutton(frame_2, text = "OK", 
                       command = ok_callback)
cancel_button <- ttkbutton(frame_2, text = "Cancel", 
                    command = function() tkdestroy(window))
#
tkpack(ok_button, side = "left", padx = 12)  # give some space
tkpack(cancel_button)


###################################################
### code chunk number 9: ex-tcltk-simple-dialog.Rnw:105-110
###################################################
tkbind("TButton", "<Return>", function(W) tcl(W, "invoke"))
tkbind("TButton", "<FocusIn>", function(W) 
       tcl(W, "state", "active"))
tkbind("TButton", "<FocusOut>", function(W) 
       tcl(W, "state", "!active"))


###################################################
### code chunk number 10: ex-tcltk-simple-dialog.Rnw:115-118
###################################################
tkwm.state(window, "normal")
tkwm.resizable(window, FALSE, FALSE)
tkfocus(ok_button)


