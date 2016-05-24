

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


###################################################
### code chunk number 201: ex-tktext-easiest
###################################################
window <- tktoplevel()
tkwm.title(window, "Simple tktext example")
txt_widget <- tktext(window)
addScrollbars(window, txt_widget)


###################################################
### code chunk number 202: tkinsert-example
###################################################
tkinsert(txt_widget, 
         "1.0", 
         paste("Lorem ipsum dolor",
               "sit amet,", sep = "\n"))


###################################################
### code chunk number 203: get-values
###################################################
value <- tkget(txt_widget, "1.0", "end")
as.character(value)                     # wrong way
tclvalue(value)


###################################################
### code chunk number 204: ScrollableWidgets.Rnw:173-174 (eval = FALSE)
###################################################
## value <- tkget(txt_widget, "1.0", "end")


###################################################
### code chunk number 205: ScrollableWidgets.Rnw:217-218
###################################################
tkinsert(txt_widget, "end", "last words", "lastWords") 


###################################################
### code chunk number 206: ScrollableWidgets.Rnw:225-227
###################################################
tktag.add(txt_widget, "first_word", 
          "1.0 wordstart", "1.0 wordend")


###################################################
### code chunk number 207: ScrollableWidgets.Rnw:231-233
###################################################
tktag.configure(txt_widget, "first_word", foreground = "red", 
                font = "helvetica 12 bold")


###################################################
### code chunk number 208: ScrollableWidgets.Rnw:250-254
###################################################
has_selection <- function(W) {
  ranges <- tclvalue(tcl(W, "tag", "ranges", "sel"))
  length(ranges) > 1 || ranges != ""
}


###################################################
### code chunk number 209: cutSelection
###################################################
tcl("tk_textCut", txt_widget)


###################################################
### code chunk number 210: ScrollableWidgets.Rnw:273-279
###################################################
popup_context <- function(W, x, y) {
  ## or use sprintf("@%s,%s", x, y) for "current"
  cur <- tkget(W, "current wordstart", "current wordend") 
  cur <- tclvalue(cur)
  popup_context_menu_for(cur, x, y)        # some function
}


###################################################
### code chunk number 211: ScrollableWidgets.Rnw:289-293
###################################################
tkmark.set(txt_widget, "leftlimit", "insert")
tkmark.gravity(txt_widget, "leftlimit", "left") # keep on left
tkinsert(txt_widget, "insert", "new text")
tkget(txt_widget, "leftlimit", "insert")


###################################################
### code chunk number 212: ScrollableWidgets.Rnw:306-308
###################################################
tcl(txt_widget, "edit", "undo")                  # no output
tcl(txt_widget, "edit", "modified")              # 1 = TRUE

