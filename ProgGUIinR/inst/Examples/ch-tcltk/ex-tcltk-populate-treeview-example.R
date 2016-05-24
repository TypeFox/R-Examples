
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
### code chunk number 240: ScrollableWidgets.Rnw:606-631
###################################################
populate_rectangular_treeview <- function(parent, m) {
  enc_frame <- ttkframe(parent)
  frame <- ttkframe(enc_frame)
  tkpack(frame, expand = TRUE, fill = "both")
  treeview <- ttktreeview(frame,
                    columns = seq_len(ncol(m)),
                    show = "headings")
  addScrollbars(frame, treeview)
  tkpack.propagate(enc_frame, FALSE)    # size from frame
  ## headings,widths
  font_measure <- tcl("font","measure","TkTextFont","0")
  charWidth <- as.integer(tclvalue(font_measure))
  sapply(seq_len(ncol(m)), function(i) {
    tcl(treeview, "heading", i, text = colnames(m)[i])
    tcl(treeview, "column", i, 
        width = 10 + charWidth*max(apply(m, 2, nchar)))
  })
  tcl(treeview, "column", ncol(m), stretch = TRUE)
  ## values
  if(ncol(m) == 1)  m <- as.matrix(paste("{", m, "}", sep=""))
  apply(m, 1, function(vals) 
    tcl(treeview, "insert", "", "end", values = vals))
  ##
  return(list(treeview = treeview, frame = enc_frame))
}


###################################################
### code chunk number 241: populate_test
###################################################
window <- tktoplevel()
m <- sapply(mtcars, as.character)
a <- populate_rectangular_treeview(window, m)
tkconfigure(a$treeview, selectmode = "extended") # multiple 
tkconfigure(a$frame, width = 300, height = 200) # frame size
tkpack(a$frame, expand = TRUE, fill = "both")
