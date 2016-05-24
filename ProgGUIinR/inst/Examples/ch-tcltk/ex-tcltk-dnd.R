### R code from vignette source 'ex-tcltk-dnd.Rnw'

###################################################
### code chunk number 1: ex-tcltk-dnd.Rnw:19-23
###################################################
## ex-tcltk-dnd
## Example of drag and drop
## idea from http://wiki.tcl.tk/416
library(tcltk)


###################################################
### code chunk number 2: ex-tcltk-dnd.Rnw:29-35
###################################################
window <- tktoplevel()
b_drag <- ttkbutton(window, text = "Drag me")
b_drop <- ttkbutton(window, text = "Drop here")
tkpack(b_drag)
tkpack(ttklabel(window, text = "Drag over me"))
tkpack(b_drop)


###################################################
### code chunk number 3: dndGlobalVariables
###################################################
.dragging <- FALSE                 # currently dragging?
.drag_value <- ""                  # value to transfer
.last_widget_id <- ""              # last widget dragged over


###################################################
### code chunk number 4: defineDragLabel
###################################################
## make a drag label to pop up during dragging to show text
## label idea from http://www.codebykevin.com/opensource/xplat_oss.html
## not shown in text, but here as an example if desired
tclServiceMode(FALSE)
tcl("image", "create", "photo", "::dnd::icon",
    data = paste(                         # uuencode gif image of document
      "R0lGODlhEAAQALMAAAAAAMbGxv//////////////////////////////////",
      "/////////////////////yH5BAEAAAEALAAAAAAQABAAAAQwMMhJ6wQ4YyuB",
      "+OBmeeDnAWNpZhWpmu0bxrKAUu57X7VNy7tOLxjIqYiapIjDbDYjADs=",
      sep = ""))
drag_label <- tktoplevel()
tkwm.title(drag_label, "")
tkwm.resizable(drag_label, FALSE,FALSE)  # no size handler in Mac OS X
tkwm.withdraw(drag_label)                # minimize window
tclServiceMode(TRUE)                    # back to redrawing
tkwm.overrideredirect(drag_label, TRUE)  # Should remove title bar (not in Mac OS X)
## can set this using tkconfigure
drag_label_label <- ttklabel(drag_label, image = "::dnd::icon", text = "this space for rent", compound = "left")
tkpack(drag_label_label)


###################################################
### code chunk number 5: buttonPressEvent
###################################################
tkbind(b_drag,"<ButtonPress-1>",function(W) {
  .dragging <<-  TRUE
  .drag_value <<- as.character(tkcget(W,text = NULL))
  .last_widget_id <<- as.character(W)
})


###################################################
### code chunk number 6: ex-tcltk-dnd.Rnw:84-86
###################################################
## If using drag_label, you can add this line to handler code:
## tkconfigure(drag_label_label, text = .drag_value)


###################################################
### code chunk number 7: ex-tcltk-dnd.Rnw:103-122
###################################################
tkbind(window, "<B1-Motion>", function(W, X, Y) {
  if(!.dragging) return()
  ## see cursor help page in API for more options
  ## For custom cursors cf. http://wiki.tcl.tk/8674. 
  tkconfigure(W, cursor = "coffee_mug")  # set cursor

  win <- tkwinfo("containing", X, Y)    # widget mouse is over
  if(as.logical(tkwinfo("exists", win))) # does widget exist?
    tkevent.generate(win, "<<DragOver>>")

  ## generate drag leave if we left last widget
  if(as.logical(tkwinfo("exists", win)) &&
     nchar(as.character(win)) > 0 && 
     length(.last_widget_id) > 0) {     # if not character(0) 
    if(as.character(win) != .last_widget_id) 
      tkevent.generate(.last_widget_id, "<<DragLeave>>")
  }
  .last_widget_id <<- as.character(win)
})


###################################################
### code chunk number 8: ex-tcltk-dnd.Rnw:125-129
###################################################
## if doing example with drag label, include the following in the motion event callback
## XXX Doesn't raise window above
##  tkwm.deiconify(drag_label)
##  tkwm.geometry(drag_label, paste("+",X,"+",Y, sep = "")) # put in offset if desired


###################################################
### code chunk number 9: ex-tcltk-dnd.Rnw:136-149
###################################################
 tkbind(b_drag, "<ButtonRelease-1>", function(W, X, Y) {
  if(!.dragging) return()
  w <- tkwinfo("containing", X, Y)
    
  if(as.logical(tkwinfo("exists", w))) {
    tkevent.generate(w, "<<DragLeave>>")
    tkevent.generate(w, "<<DragDrop>>")
    tkconfigure(w, cursor = "")
  }
  .dragging <<- FALSE
  .last_widget_id <<- "" 
  tkconfigure(W, cursor = "")
})


###################################################
### code chunk number 10: ex-tcltk-dnd.Rnw:152-154
###################################################
## if doing drag_label, we have this as well:
## tkwm.withdraw(drag_label)


###################################################
### code chunk number 11: ex-tcltk-dnd.Rnw:161-165
###################################################
tkbind(b_drop,"<<DragOver>>",function(W) {
  if(.dragging) 
    tcl(W, "state", "active")
})


###################################################
### code chunk number 12: ex-tcltk-dnd.Rnw:170-176
###################################################
tkbind(b_drop,"<<DragLeave>>", function(W) {
  if(.dragging)  {
    tkconfigure(W, cursor = "")
    tcl(W, "state", "!active")  
   }
})


###################################################
### code chunk number 13: ex-tcltk-dnd.Rnw:182-186
###################################################
tkbind(b_drop,"<<DragDrop>>", function(W) {
  tkconfigure(W, text = .drag_value)
  .drag_value <- ""
})


