
###################################################
### code chunk number 63: newWindow
###################################################
newWindow <- function(title, command, parent,
                      width, height) {
  window <- tktoplevel()

  if(!missing(title)) tkwm.title(window, title)

  if(!missing(command)) 
    tkwm.protocol(window, "WM_DELETE_WINDOW", function() {
      if(command())            # command returns logical
        tkdestroy(window)
    })

  if(!missing(parent)) {
    parent_window <- tkwinfo("top-level", parent)
    if(as.logical(tkwinfo("viewable", parent_window))) {
      tkwm.transient(window, parent)
    }
  }
  
  if(!missing(width)) tkconfigure(window, width = width)
  if(!missing(height)) tkconfigure(window, height = height)

  window_system <- tclvalue(tcl("tk", "windowingsystem"))
  if(window_system == "aqua") {
    frame <- ttkframe(window, padding = c(3,3,12,12))
  } else {
    int_frame <- ttkframe(window, padding = 0)
    tkpack(int_frame, expand = TRUE, fill = "both")
    frame <- ttkframe(int_frame, padding = c(3,3,12,0))
    sizegrip <- ttksizegrip(int_frame)
    tkpack(sizegrip, side = "bottom", anchor = "se")
  }
  tkpack(frame, expand = TRUE, fill = "both", side = "top")

  return(frame)
}
