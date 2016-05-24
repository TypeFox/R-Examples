
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
### code chunk number 235: treeExample
###################################################
window <- tktoplevel()
tkwm.title(window, "Choose CRAN mirror")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 236: ScrollableWidgets.Rnw:497-503
###################################################
treeview <- 
  ttktreeview(frame, 
              columns = 1,        # column identifier is "1"
              show = "headings",  # not "#0"
              height = 25)        
addScrollbars(frame, treeview)    # our scroll bar function


###################################################
### code chunk number 237: getCRANmirrors
###################################################
x <- getCRANmirrors()
Host <- x$Host
shade <- c("none", "gray")                     # tag names
for(i in seq_along(Host))
  ID <- tkinsert(treeview, "", "end", 
                 values = as.tclObj(Host[i]),
                 tag = shade[i %% 2])          # none or gray
tktag.configure(treeview, "gray", background = "gray95") 


###################################################
### code chunk number 238: ScrollableWidgets.Rnw:579-580
###################################################
tcl(treeview, "heading", 1, text = "Host", anchor = "center")


###################################################
### code chunk number 239: ScrollableWidgets.Rnw:595-597
###################################################
tcl(treeview, "column", 1, width = 400,  
    stretch = TRUE, anchor = "w")
