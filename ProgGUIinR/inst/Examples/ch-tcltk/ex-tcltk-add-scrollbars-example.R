
###################################################
### code chunk number 200: addScrollbars
###################################################
## function to add scroll bars to widget and pack into grid
## parent uses grid manager -- con't pack in other children
addScrollbars <- function(parent, widget) {
  xscr <- ttkscrollbar(parent, orient = "horizontal",
                       command = function(...) tkxview(widget, ...))
  yscr <- ttkscrollbar(parent, orient = "vertical",
                       command = function(...) tkyview(widget, ...))
  ##
  tkconfigure(widget,
              xscrollcommand = function(...) tkset(xscr,...),
              yscrollcommand = function(...) tkset(yscr,...))
  ##
  tkgrid(widget, row = 0, column = 0, sticky = "news")
  tkgrid(yscr, row = 0,column = 1, sticky = "ns")
  tkgrid(xscr, row = 1, column = 0, sticky = "ew")
  tkgrid.columnconfigure(parent, 0, weight = 1)
  tkgrid.rowconfigure(parent, 0, weight = 1)
}

