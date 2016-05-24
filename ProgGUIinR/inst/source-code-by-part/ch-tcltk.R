### R code from vignette source 'ch-tcltk.Rnw'

###################################################
### code chunk number 1: ch-tcltk.Rnw:3-28
###################################################
require(tcltk)
source("../booktabs.R")

noPrompt <- function() {
  options(prompt=" ")
  options(continue=" ")
}
doPrompt <- function() {
  options(prompt="> ")
  options(continue="+ ")
}
noPrompt()
noTclPrint <- function(x, ...) {
  y <- try(as.character(x))
  if(inherits(y,"try-error")) {
    cat(y)
  } else if(length(y) == 0 ||
            (length(y) == 1 && nchar(y) == 0 )) {
    return()
    invisible(x)
  }
  tcltk:::print.tclObj(x)
}
print.tclObj <- noTclPrint



###################################################
### code chunk number 2: simpleExample
###################################################
library(tcltk)
##
window <- tktoplevel()
tkwm.title(window, "Simple dialog")
##
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")
##
nested_frame <- ttkframe(frame); tkpack(nested_frame)
##
label <- ttklabel(nested_frame, text = "Enter your name:")
tkpack(label, side = "left")
##
text_var <- tclVar("")
entry <- ttkentry(nested_frame, textvariable = text_var)
tkpack(entry)
##
button_frame <- ttkframe(frame)
tkpack(button_frame, anchor = "ne")
button <- ttkbutton(button_frame, text = "Click")
tkpack(button, side = "right")
##
handler <- function() {
  msg <- sprintf("Hello %s", tclvalue(text_var))
  print(msg)
}
tkconfigure(button, command = handler)


###################################################
### code chunk number 3: Overview.Rnw:170-173
###################################################
library(tcltk)
.Tcl("set x {some text}")               # assignment
.Tcl("puts $x")                         # prints to stdout


###################################################
### code chunk number 4: Overview.Rnw:181-182
###################################################
.Tcl("string length $x")                # call a command


###################################################
### code chunk number 5: Overview.Rnw:210-212
###################################################
## To simplify coercion to logical, we define a new method:
as.logical.tclObj <- function(x, ...) as.logical(as.numeric(x))


###################################################
### code chunk number 6: Overview.Rnw:227-231
###################################################
## set up label for Sweave below. Not shown
window <- tktoplevel()
label <- ttklabel(window)
tkpack(label)


###################################################
### code chunk number 7: Overview.Rnw:234-235
###################################################
tkconfigure(label, text = "new text")


###################################################
### code chunk number 8: eval
###################################################
window <- tktoplevel()


###################################################
### code chunk number 9: Overview.Rnw:350-351
###################################################
label <- ttklabel(nested_frame, text = "Enter your name:")


###################################################
### code chunk number 10: Overview.Rnw:389-390
###################################################
head(as.character(tkconfigure(label)))      # first 6 only


###################################################
### code chunk number 11: Overview.Rnw:400-402
###################################################
tkcget(label, "-text")               # retrieve text property
tkcget(label, text = NULL)           # alternate syntax


###################################################
### code chunk number 12: Overview.Rnw:409-410
###################################################
as.character(tkcget(label, text = NULL))


###################################################
### code chunk number 13: Overview.Rnw:414-415
###################################################
tclvalue(tkcget(label, text = NULL))


###################################################
### code chunk number 14: str
###################################################
str(button)


###################################################
### code chunk number 15: Overview.Rnw:453-457 (eval = FALSE)
###################################################
## tkpack(frame, expand = TRUE, fill = "both")
## tkpack(label, side = "left")
## tkpack(entry)
## tkpack(button_frame, anchor = "ne")


###################################################
### code chunk number 16: Overview.Rnw:485-487 (eval = FALSE)
###################################################
## text_var <- tclVar("")
## entry <- ttkentry(g, textvariable = text_var)


###################################################
### code chunk number 17: Overview.Rnw:502-504
###################################################
tclvalue(text_var) <- "Somebody's name"
tclvalue(text_var)


###################################################
### code chunk number 18: Overview.Rnw:508-509
###################################################
as.character(text_var)


###################################################
### code chunk number 19: tclVarExample (eval = FALSE)
###################################################
## x <- tclVar(1)
## print(x)                           # a tclVar object
## tclvalue(x)                        # internal name
## as.numeric(tclvalue(x))            # get value then coerce


###################################################
### code chunk number 20: tclArrayExample
###################################################
x <- tclArray()                    # no init
x$one <- 1; x[[2]] <- 2            # $<- and [[<-
x[[1]]                             # no match by index
names(x)                           # the stored names
x[['2']]                           # match by name, not index


###################################################
### code chunk number 21: Overview.Rnw:550-557
###################################################
button <- ttkbutton(button_frame, text = "Click")
#
handler <- function() {
  msg <- sprintf("Hello %s", tclvalue(text_var))
  print(msg)
}
tkconfigure(button, command = handler)


###################################################
### code chunk number 22: Overview.Rnw:594-597 (eval = FALSE)
###################################################
## window <- tktoplevel()
## frame <- ttkframe(window, padding = c(3,3,12,12))
## tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 23: Overview.Rnw:626-627
###################################################
.Tcl("ttk::style theme names")


###################################################
### code chunk number 24: Overview.Rnw:631-632
###################################################
.Tcl("ttk::style theme use clam")


###################################################
### code chunk number 25: has-focus
###################################################
as.logical(tcl(button, "instate", "focus"))


###################################################
### code chunk number 26: Overview.Rnw:646-647
###################################################
tcl(button, "state", "disabled")             # not sensitive


###################################################
### code chunk number 27: Overview.Rnw:650-651
###################################################
tcl(button, "state", "!disabled")            # sensitive again


###################################################
### code chunk number 28: WindowClass
###################################################
tkwinfo("class", label)


###################################################
### code chunk number 29: Overview.Rnw:756-759
###################################################
window <- tktoplevel()
ttkbutton(window)
ttkbutton(window)


###################################################
### code chunk number 30: classOfChildren
###################################################
(children <- tkwinfo("children", window))
sapply(as.character(children), function(i) tkwinfo("class",i))


###################################################
### code chunk number 31: Overview.Rnw:796-798
###################################################
tkconfigure(label, foreground = "red")
tkconfigure(label, foreground = "#00aa00")


###################################################
### code chunk number 32: Overview.Rnw:819-820
###################################################
tkconfigure(label, font = "TkFixedFont")


###################################################
### code chunk number 33: Overview.Rnw:827-841
###################################################
DF <- rbind(
            c("\\code{TkDefaultFont}","Default font for all GUI items not otherwise specified"),
            c("\\code{TkTextFont}","Font for text widgets"),
            c("\\code{TkFixedFont}","Fixed-width font"),
            c("\\code{TkMenuFont}","Menu bar fonts"),
            c("\\code{TkHeadingFont}","Font for column headings"),
            c("\\code{TkCaptionFont}","Caption font (dialogs)"),
            c("\\code{TkSmallCaptionFont}","Smaller caption font"),
            c("\\code{TkIconFont}","Icon and text font")
            )
colnames(DF) <- c("Standard font name", "Description")
cat(booktabs(DF, 
             caption="Standard font names defined by a theme.",
             label="tab:tcltk-std-fonts"))


###################################################
### code chunk number 34: Overview.Rnw:846-849
###################################################
tkfont.create("our_font", family = "Helvetica", size = 12, 
              weight = "bold")
tkconfigure(label, font = "our_font")


###################################################
### code chunk number 35: fontSizes
###################################################
font_measure <- tcl("font", "measure", "TkTextFont", "M")
font_width <- as.integer(tclvalue(font_measure))
tmp <- tkfont.metrics("TkTextFont", "linespace" = NULL)
font_height <- as.numeric(tclvalue(tmp))
#
c(width = font_width, height = font_height)


###################################################
### code chunk number 36: Overview.Rnw:912-913
###################################################
tkimage.create("photo", "::img::tclLogo", file = "tclp.gif")


###################################################
### code chunk number 37: Overview.Rnw:938-940
###################################################
label <- ttklabel(window, image = "::img::tclLogo", 
                  text = "logo text", compound = "top")


###################################################
### code chunk number 38: Overview.Rnw:942-943
###################################################
tkpack(label)


###################################################
### code chunk number 39: Overview.Rnw:958-959
###################################################
tkconfigure("::img::tclLogo", palette = 16)


###################################################
### code chunk number 40: Overview.Rnw:964-966
###################################################
## redo changing of palette
tkconfigure("::img::tclLogo", palette = "fullcolor")


###################################################
### code chunk number 41: Ctrl-q-binding
###################################################
window <- tktoplevel()
button <- ttkbutton(window, text = "Some widget with focus")
tkpack(button)
tkbind(window, "<Control-q>", function() tkdestroy(window))


###################################################
### code chunk number 42: bind_examples
###################################################
window <- tktoplevel()
label <- ttklabel(window, text = "Click Ok for a message")
button1 <- ttkbutton(window, text = "Cancel", 
                command = function() tkdestroy(window))
button2 <- ttkbutton(window, text = "Ok", command=function() {
  print("initiate an action")
})
sapply(list(label, button1, button2), tkpack)
##
tkbind(window, "C", function() tcl(button1, "invoke"))
tkconfigure(button1, underline = 0)
##
tkbind(window, "O", function() tcl(button1, "invoke"))
tkconfigure(button2, underline = 0)
tkfocus(button2)
##
tkbind("TButton", "<Return>", function(W) {
  tcl(W, "invoke")
})


###################################################
### code chunk number 43: Overview.Rnw:1112-1114
###################################################
 tkevent.add("<<Paste>>", "<Control-y>")
 tkevent.add("<<Save>>", "<Control-x><Control-s>")


###################################################
### code chunk number 44: Overview.Rnw:1123-1124
###################################################
tkevent.generate(button, "<<Save>>")


###################################################
### code chunk number 45: Overview.Rnw:1192-1201
###################################################
window <- tktoplevel()
button <- 
  ttkbutton(window, text = "Click me for the x,y position")
tkpack(button)
tkbind(button, "<ButtonPress-1>", function(W, x, y, X, Y) {
  print(W)                              # an ID
  print(c(x, X))                        # character class
  print(c(y, Y))
  })


###################################################
### code chunk number 46: ReturnValueDoesNotStopPropogation (eval = FALSE)
###################################################
## ## This example illustrates how handlers are different in tcltk than Tk
## ## In this case we can't break out of the chain of handlers called in tcltk
## library(tcltk)
## window <- tktoplevel()
## button <- ttkbutton(window, text = "click me", command = function() {
##   ## default command: bind Button <1> {tk::ButtonDown %W}
##   print("hi")
## })
## 
## tkpack(button)
## 
## tkbind(button, "<Button-1>", function() {
##   print("hello")
##   .Tcl("return -code break") ## doesn't suppress call of other
## })
## 
## ## in Tcl returning -code break suppresses call to next handler
## .Tcl("toplevel .w")
## .Tcl("wm title .w {title}")
## .Tcl("button .w.b -text {Click me} -command {puts hi}")
## .Tcl("pack .w.b")
## #.Tcl("bind .w.b <Button-1> {puts hello}")
## .Tcl("bind .w.b <Button-1> {puts hello; return -code break}") ## no hi


###################################################
### code chunk number 47: Overview.Rnw:1246-1261 (eval = FALSE)
###################################################
## after_ID <- ""
## some_flag <- TRUE
## repeat_call <- function(ms = 100, f) {
##   after_ID <<- tcl("after", ms, function() {
##     if(someFlag) {                      
##       f()
##       after_ID <<- repeat_call(ms, f)
##     }  else {
##       tcl("after", "cancel", after_ID)
##     }
##   })
## }
## repeat_call(2000, function() {
##   print("Running. Set someFlag <- FALSE to stop.")
## })


###################################################
### code chunk number 48: ex-tcltk-dnd.Rnw:19-23
###################################################
## ex-tcltk-dnd
## Example of drag and drop
## idea from http://wiki.tcl.tk/416
library(tcltk)


###################################################
### code chunk number 49: ex-tcltk-dnd.Rnw:29-35
###################################################
window <- tktoplevel()
b_drag <- ttkbutton(window, text = "Drag me")
b_drop <- ttkbutton(window, text = "Drop here")
tkpack(b_drag)
tkpack(ttklabel(window, text = "Drag over me"))
tkpack(b_drop)


###################################################
### code chunk number 50: dndGlobalVariables
###################################################
.dragging <- FALSE                 # currently dragging?
.drag_value <- ""                  # value to transfer
.last_widget_id <- ""              # last widget dragged over


###################################################
### code chunk number 51: defineDragLabel
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
### code chunk number 52: buttonPressEvent
###################################################
tkbind(b_drag,"<ButtonPress-1>",function(W) {
  .dragging <<-  TRUE
  .drag_value <<- as.character(tkcget(W,text = NULL))
  .last_widget_id <<- as.character(W)
})


###################################################
### code chunk number 53: ex-tcltk-dnd.Rnw:84-86
###################################################
## If using drag_label, you can add this line to handler code:
## tkconfigure(drag_label_label, text = .drag_value)


###################################################
### code chunk number 54: ex-tcltk-dnd.Rnw:103-122
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
### code chunk number 55: ex-tcltk-dnd.Rnw:125-129
###################################################
## if doing example with drag label, include the following in the motion event callback
## XXX Doesn't raise window above
##  tkwm.deiconify(drag_label)
##  tkwm.geometry(drag_label, paste("+",X,"+",Y, sep = "")) # put in offset if desired


###################################################
### code chunk number 56: ex-tcltk-dnd.Rnw:136-149
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
### code chunk number 57: ex-tcltk-dnd.Rnw:152-154
###################################################
## if doing drag_label, we have this as well:
## tkwm.withdraw(drag_label)


###################################################
### code chunk number 58: ex-tcltk-dnd.Rnw:161-165
###################################################
tkbind(b_drop,"<<DragOver>>",function(W) {
  if(.dragging) 
    tcl(W, "state", "active")
})


###################################################
### code chunk number 59: ex-tcltk-dnd.Rnw:170-176
###################################################
tkbind(b_drop,"<<DragLeave>>", function(W) {
  if(.dragging)  {
    tkconfigure(W, cursor = "")
    tcl(W, "state", "!active")  
   }
})


###################################################
### code chunk number 60: ex-tcltk-dnd.Rnw:182-186
###################################################
tkbind(b_drop,"<<DragDrop>>", function(W) {
  tkconfigure(W, text = .drag_value)
  .drag_value <- ""
})


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


###################################################
### code chunk number 64: rotateLabels
###################################################
window <- tktoplevel()
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")
#
x <- strsplit("Lorem ipsum dolor sit amet ...", "\\s")[[1]]
labels <- lapply(x, function(i) ttklabel(frame, text = i))
sapply(labels, function(i) tkpack(i, side = "left"))
#
rotateLabel <- function() {
  children <- as.character(tkpack("slaves", frame))
  tkpack.forget(children[1])
  tkpack(children[1], after = children[length(children)], 
         side = "left")
}


###################################################
### code chunk number 65: BasicContainers.Rnw:341-342 (eval = FALSE)
###################################################
## for(i in 1:20) {rotateLabel(); Sys.sleep(1)}


###################################################
### code chunk number 66: padding-pady-ipady
###################################################
## Code to make padding/pady/ipady figure
## padding
w <- tktoplevel();
tkwm.title(w, "Padding/pady/ipady")
f1 <- ttkframe(w, padding = c(3,3,12,12))  # Some breathing room
tkpack(f1, expand = TRUE, fill = "both")     # For window expansion


f <- ttkframe(f1,padding = 20, border = 1);
tkpack(f, side = "left", expand = TRUE, fill = "x")
tkpack(ttkbutton(f, text = "padding"))
tkpack(ttkbutton(f, text = "padding"))

tkpack(ttkseparator(f1, orient = "vertical"), expand = TRUE, fill = "y", side = "left")

##  pady outside border
f <- ttkframe(f1, border = 1)
tkpack(f, side = "left", expand = TRUE, fill = "both")
tkpack(ttkbutton(f, text = "pady"), pady = 10)
tkpack(ttkbutton(f, text = "pady"), pady = 10)

tkpack(ttkseparator(f1, orient = "vertical"), expand = TRUE, fill = "y", side = "left")
## ipady within border
f <- ttkframe(f1, border = 1)
tkpack(f, side = "left", expand = TRUE, fill = "both")
tkpack(ttkbutton(f, text = "ipady"), ipady = 10)
tkpack(ttkbutton(f, text = "ipady"), ipady = 10)


###################################################
### code chunk number 67: BasicContainers.Rnw:404-416
###################################################
## not shown
## sets up window with some child components to illustrate
## tkpack("configure", ...)
window <- tktoplevel()
tkwm.title(window, "tkpack configure example")
frame <- ttkframe(window, padding = 5)
tkpack(frame, expand = TRUE, fill = "both")
l <- list()
for(i in 1:5) {
  l[[i]] <- ttklabel(frame, text = i)
  tkpack(l[[i]], side = "left")
}


###################################################
### code chunk number 68: BasicContainers.Rnw:423-426
###################################################
all_children <- as.character(tkwinfo("children", frame))
sapply(all_children, tkpack.configure,  padx = 10, pady = 5)



###################################################
### code chunk number 69: anchorPoints
###################################################
## Example of anchor arguments
## not shown in text
w <- tktoplevel()
tkwm.title(w, "anchor arguments")
f <- ttkframe(w, padding = c(3,3,12,12)); tkpack(f, expand = TRUE, fill = "both")
tr <- ttkframe(f); mr <- ttkframe(f); br <- ttkframe(f)
sapply(list(tr, mr, br), tkpack, expand = TRUE, fill = "both")
tkpack(ttklabel(tr, text = "nw"), anchor = "nw", expand = TRUE, side = "left")
tkpack(ttklabel(tr, text = "n"),  anchor = "n",  expand = TRUE, side = "left")
tkpack(ttklabel(tr, text = "ne"), anchor = "ne", expand = TRUE, side = "left")
##
tkpack(ttklabel(mr, text = "w"),       anchor = "w",       expand = TRUE, side = "left")
tkpack(ttklabel(mr, text = "center"),  anchor = "center",  expand = TRUE, side = "left")
tkpack(ttklabel(mr, text = "e"),       anchor = "e",       expand = TRUE, side = "left")
##
tkpack(ttklabel(br, text = "sw"), anchor = "sw", expand = TRUE,  side = "left")
tkpack(ttklabel(br, text = "s"),  anchor = "s",  expand = TRUE,  side = "left")
tkpack(ttklabel(br, text = "se"), anchor = "se", expand = TRUE,  side = "left")


###################################################
### code chunk number 70: expandFill
###################################################
window <- tktoplevel()
tkwm.title(window, "Expand/Fill arguments")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")
##
pack_btn <- function(txt, ...) 
  tkpack(button <- ttkbutton(frame, text = txt), ...)
##
pack_btn("Top",    side="top",    expand=TRUE, fill="both") 
pack_btn("Bottom", side="bottom", expand=TRUE, fill="both") 
pack_btn("Left",   side="left",   expand=TRUE, fill="both") 
pack_btn("Right",  side="right",  expand=TRUE, fill="both") 


###################################################
### code chunk number 71: BasicContainers.Rnw:537-539
###################################################
children <- as.character(tkwinfo("children", frame))
sapply(children, tkpack.configure, fill = "none")


###################################################
### code chunk number 72: BasicContainers.Rnw:554-555 (eval = FALSE)
###################################################
## tkwm.geometry(window, "")


###################################################
### code chunk number 73: BasicContainers.Rnw:560-561 (eval = FALSE)
###################################################
## tkwinfo("top-level", button)


###################################################
### code chunk number 74: NotEnoughSpace
###################################################
## example of last in first covered
## Create this GUI, then shrink window with the mouse
w <- tktoplevel()
f <- ttkframe(w); tkpack(f, expand = TRUE, fill = "both")

g1 <- ttkframe(f); tkpack(g1, expand = TRUE, fill = "both")
g2 <- ttkframe(f); tkpack(g2, expand = TRUE, fill = "both")

b11 <- ttkbutton(g1, text = "first")
b12 <- ttkbutton(g1, text = "second")
b21 <- ttkbutton(g2, text = "first")
b22 <- ttkbutton(g2, text = "second")

tkpack(b11, side = "left"); tkpack(b12, side = "left")
tkpack(b21, side="right"); tkpack(b22, side = "right")
## Now shrink window



###################################################
### code chunk number 75: ex-tcltk-pack.Rnw:1-4
###################################################
library(tcltk)
## pack examples
## how to pack in buttons


###################################################
### code chunk number 76: ex-tcltk-pack.Rnw:7-11
###################################################
window <- tktoplevel()
tkwm.title(window, "Examples using pack as a geometry manager")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 77: ex-tcltk-pack.Rnw:31-36
###################################################
frame_1 <- ttklabelframe(frame, text="plain vanilla")
tkpack(frame_1, expand = TRUE, fill = "x")
l <- function(f) 
  list(ttkbutton(f, text="cancel"), ttkbutton(f, text="ok"))
sapply(l(frame_1), tkpack, side = "left")


###################################################
### code chunk number 78: moveRight
###################################################
frame_2 <- ttklabelframe(frame, text = "push to right")
tkpack(frame_2, expand = TRUE, fill = "x")
tkpack(ttklabel(frame_2, text = " "), 
       expand = TRUE, fill = "x", side = "left")
sapply(l(frame_2), tkpack, side = "left")


###################################################
### code chunk number 79: appleSays
###################################################
frame_3 <- ttklabelframe(frame, text="push right with space")
tkpack(frame_3, expand = TRUE, fill = "x")
tkpack(ttklabel(frame_3, text = " "), expand=TRUE, fill="x", 
       side = "left")
sapply(l(frame_3), tkpack, side = "left", padx = 6)


###################################################
### code chunk number 80: ex-tcltk-simple-dialog.Rnw:11-15
###################################################
## example of  simple non-modal dialog
## much taken from msgbox.tcl in tk source

library(tcltk)


###################################################
### code chunk number 81: ex-tcltk-simple-dialog.Rnw:19-25
###################################################
title <- "message dialog"
message <- "Do you like tcltk so far?"
parent <- NULL
tkimage.create("photo", "::img::tclLogo", 
               file = system.file("images","tclp.gif",
                 package = "ProgGUIinR"))


###################################################
### code chunk number 82: ex-tcltk-simple-dialog.Rnw:30-35
###################################################
window <- tktoplevel()
tkwm.title(window, title)
tkwm.state(window, "withdrawn")
frame <- ttkframe(window,  padding = c(3, 3, 12, 12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 83: ex-tcltk-simple-dialog.Rnw:45-58
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
### code chunk number 84: ex-tcltk-simple-dialog.Rnw:63-65
###################################################
label <- ttklabel(frame, image = "::img::tclLogo", padding=5) 
tkpack(label, side = "left")


###################################################
### code chunk number 85: ex-tcltk-simple-dialog.Rnw:71-76
###################################################
frame_1 <- ttkframe(frame)
tkpack(frame_1, expand = TRUE, fill = "both")
#
m <- ttklabel(frame_1, text = message)
tkpack(m, expand = TRUE, fill = "both")


###################################################
### code chunk number 86: ex-tcltk-simple-dialog.Rnw:80-82
###################################################
frame_2 <- ttkframe(frame_1)
tkpack(frame_2)


###################################################
### code chunk number 87: ex-tcltk-simple-dialog.Rnw:86-97
###################################################
ok_callback <- function() {
  print("That's great")
  tkdestroy(window)
}
ok_button <- ttkbutton(frame_2, text = "OK", 
                       command = ok_callback)
cancel_button <- ttkbutton(frame_2, text = "Cancel", 
                    command = function() tkdestroy(window))
#
tkpack(ok_button, side = "left", padx = 12)  # give some space
tkpack(cancel_button)


###################################################
### code chunk number 88: ex-tcltk-simple-dialog.Rnw:105-110
###################################################
tkbind("TButton", "<Return>", function(W) tcl(W, "invoke"))
tkbind("TButton", "<FocusIn>", function(W) 
       tcl(W, "state", "active"))
tkbind("TButton", "<FocusOut>", function(W) 
       tcl(W, "state", "!active"))


###################################################
### code chunk number 89: ex-tcltk-simple-dialog.Rnw:115-118
###################################################
tkwm.state(window, "normal")
tkwm.resizable(window, FALSE, FALSE)
tkfocus(ok_button)


###################################################
### code chunk number 90: BasicContainers.Rnw:639-644
###################################################
window <- tktoplevel()
tkwm.title(window, "tkgrid.rowconfigure example")
parent <- ttkframe(window); tkpack(parent)
child <- ttklabel(parent, text = "test")
tkgrid(child)


###################################################
### code chunk number 91: BasicContainers.Rnw:647-648
###################################################
tkgrid.rowconfigure(parent, 0, weight = 1)


###################################################
### code chunk number 92: ex-tcltk-toolbar.Rnw:1-4
###################################################
## Toolbar example in Tk
## tbFrame and Frame are key containers
library(tcltk)


###################################################
### code chunk number 93: ex-tcltk-toolbar.Rnw:18-21
###################################################
window <- tktoplevel(); tkwm.title(window, "Toolbar example")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 94: ex-tcltk-toolbar.Rnw:27-29
###################################################
tool_bar_frame <- ttkframe(frame, padding = 0)
content_frame <- ttkframe(frame, padding = 4)


###################################################
### code chunk number 95: ex-tcltk-toolbar.Rnw:36-44
###################################################
tkgrid(tool_bar_frame, row = 0, column = 0, sticky = "we")
tkgrid(content_frame, row = 1, column = 0, sticky  =  "news")
tkgrid.rowconfigure(frame, 0, weight = 0)
tkgrid.rowconfigure(frame, 1, weight = 1)
tkgrid.columnconfigure(frame, 0, weight = 1)
#
txt <- "Lorem ipsum dolor sit amet..." # sample text
tkpack(ttklabel(content_frame, text = txt))


###################################################
### code chunk number 96: ex-tcltk-toolbar.Rnw:50-52
###################################################
tcl("ttk::style", "configure", "Toolbar.TButton", 
    font = "helvetica 12", padding = 0, borderwidth = 0)


###################################################
### code chunk number 97: ex-tcltk-toolbar.Rnw:58-75
###################################################
make_icon <- function(parent, stock_name, command = NULL) {
  icon_file <- system.file("images", 
                          paste(stock_name,"gif",sep = "."), 
                          package = "gWidgets")
  if(nchar(icon_file) == 0) {
    b <- ttkbutton(parent, text = stock_name, width = 6)
  } else {
    icon_name <- paste("::img::",stock_name, sep = "")
    tkimage.create("photo", icon_name, file = icon_file)
    b <- ttkbutton(parent, image = icon_name, 
                   style = "Toolbar.TButton", text=stock_name, 
                   compound = "top", width = 6)
    if(!is.null(command))
      tkconfigure(b, command = command)
  }
  return(b)
}


###################################################
### code chunk number 98: addButtons
###################################################
sapply(c("ok", "quit", "cancel"), function(i)
       tkpack(make_icon(tool_bar_frame, i), side = "left"))


###################################################
### code chunk number 99: ex-tcltk-toolbar.Rnw:90-93
###################################################
setState <- function(W, state) tcl(W, "state", state)
tkbind("TButton","<Enter>",function(W) setState(W, "active"))
tkbind("TButton","<Leave>",function(W) setState(W, "!active"))


###################################################
### code chunk number 100: checkStyle
###################################################
function(W) {
  if(as.character(tkcget(W, "-style")) == "Toolbar.TButton")
    cat("... do something for toolbar buttons ...")
}


###################################################
### code chunk number 101: ex-tcltk-grid-layout.Rnw:7-34
###################################################
library(tcltk)
##' from chron with slight change to arguments
day.of.week <- function (year, month, day) {
    ix <- year + trunc((month - 14)/12)
    jx <- (trunc((13 * (month + 10 - (month + 10)%/%13 * 12) - 
        1)/5) + day + 77 + (5 * (ix - (ix%/%100) * 100))%/%4 + 
        ix%/%400 - (ix%/%100) * 2)
    jx%%7
}


## is this a valid date
validDate <- function(year, month, day) 
  !is.na(as.Date(sprintf("%s-%s-%s", year, month, day), format = "%Y-%m-%d"))

## how many days in a month
days.in.month <- function(year, month) {
  for(i in c(31, 30, 29, 28)) {
    if(validDate(year, month, i))
      return(i)
  }
}
## 0-based week of month
week.of.month <- function(year, month, day) {
  first.day <- day.of.week(year, month, 1)
  (first.day + day - 1) %/% 7
}


###################################################
### code chunk number 102: ex-tcltk-grid-layout.Rnw:47-64
###################################################
make_month <- function(parent, year, month) {
  ## add headers
  days <- c("S","M","T","W","Th","F","S")
  sapply(1:7, function(i) {
    label <- ttklabel(parent, text = days[i])       
    tkgrid(label, row = 0, column = i-1, sticky = "")
  })
  ## add days
  sapply(seq_len(ProgGUIinR:::days.in.month(year, month)), 
         function(day) {
           label <- ttklabel(parent, text = day)
           row <- ProgGUIinR:::week.of.month(year, month, day)
           col <- ProgGUIinR:::day.of.week(year, month, day)
           tkgrid(label, row = 1 + row, column = col,
                  sticky = "e")
         })
}


###################################################
### code chunk number 103: ex-tcltk-grid-layout.Rnw:69-70
###################################################
year <- 2000; month <- 1


###################################################
### code chunk number 104: ex-tcltk-grid-layout.Rnw:75-82
###################################################
window <- tktoplevel()
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both", side = "top")
c_frame <- ttkframe(frame)
cal_frame <- ttkframe(frame)
tkpack(c_frame, fill = "x", side = "top")
tkpack(cal_frame, expand = TRUE, anchor = "n")


###################################################
### code chunk number 105: ex-tcltk-grid-layout.Rnw:88-96
###################################################
previous_button <- ttklabel(c_frame, text = "<")
next_button <- ttklabel(c_frame, text = ">")
current_month <- ttklabel(c_frame)
#
tkpack(previous_button, side = "left", anchor = "w")
tkpack(current_month, side = "left", anchor = "center", 
       expand = TRUE)
tkpack(next_button, side = "left", anchor = "e")


###################################################
### code chunk number 106: stackedWidget
###################################################
set_month <- function() {
  tkpack("forget", cal_frame)
  cal_frame <<- ttkframe(frame)
  make_month(cal_frame, year, month)
  tkconfigure(current_month,              # month label
              text = sprintf("%s %s", month.abb[month], year))
  tkpack(cal_frame)
}
set_month()                              # initial calendar


###################################################
### code chunk number 107: connectSignal
###################################################
tkbind(previous_button, "<Button-1>", function() {
  if(month > 1) {
    month <<- month - 1
  } else {
    month <<- 12; year <<- year - 1
  }
  set_month()
})


###################################################
### code chunk number 108: ex-tcltk-grid-layout.Rnw:132-140
###################################################
tkbind(next_button, "<Button-1>", function() {
 if(month < 12) {
    month <<- month + 1
  } else {
    month <<- 1; year <<- year + 1
  }
  set_month()
})


###################################################
### code chunk number 109: ex-tcltk-grid-layout.Rnw:146-151
###################################################
tkbind("TLabel", "<Button-1>", function(W) {
  day <- as.numeric(tkcget(W, "-text"))
  if(!is.na(day))
    print(sprintf("You selected: %s/%s/%s", month, day, year))
})


###################################################
### code chunk number 110: panedWindowExample
###################################################
window <- tktoplevel()
tkwm.title(window, "Paned window example")
paned <- ttkpanedwindow(window, orient = "horizontal")
tkpack(paned, expand = TRUE, fill = "both")
left <- ttklabel(paned, text = "left")
right <- ttklabel(paned, text = "right")
#
tkadd(paned, left, weight = 1)
tkadd(paned, right, weight = 2)


###################################################
### code chunk number 111: BasicContainers.Rnw:732-734
###################################################
tcl(paned, "forget", right)
tkadd(paned, right, weight = 2) ## how to add back


###################################################
### code chunk number 112: BasicContainers.Rnw:743-745
###################################################
width <- as.integer(tkwinfo("width", paned))  # or "height"
tcl(paned, "sashpos", 0, floor(0.75*width))


###################################################
### code chunk number 113: BasicContainers.Rnw:754-758
###################################################
## ttknotebook example
window <- tktoplevel();
tkwm.title(window, "Notebook example")
notebook <- ttknotebook(window); tkpack(notebook, expand = TRUE, fill = "both")


###################################################
### code chunk number 114: notebookExample
###################################################
icon_file <- system.file("images",paste("help","gif",sep="."),
                        package = "gWidgets")
icon_name <- "::tcl::helpIcon"
tkimage.create("photo", icon_name, file = icon_file)
#
page2_label <- ttklabel(notebook, text = "Page 2")
tkadd(notebook, page2_label, sticky = "nswe", text="label 2", 
    image = icon_name, compound = "right")
## put page 1 label first (a tabID of 0); use tkinsert
page1_label <- ttklabel(notebook, text = "Page 1")
tkinsert(notebook, 0, page1_label, sticky = "nswe", 
         text = "label 1")


###################################################
### code chunk number 115: BasicContainers.Rnw:817-822
###################################################
tcl(notebook, "index", "current")    # current page for tabID
length(as.character(tcl(notebook,"tabs")))  # number of pages
tcl(notebook, "select", 0)           # select by index
tcl(notebook, "forget", page1_label) # "forget" removes a page
tcl(notebook, "add", page1_label)    # can be managed again.


###################################################
### code chunk number 116: notebookTraversal
###################################################
tcl("ttk::notebook::enableTraversal", notebook)


###################################################
### code chunk number 117: tkMessageBox (eval = FALSE)
###################################################
## tkmessageBox(title = "Confirm", message = "Really exit?", 
##              detail = "Changes need saving.", 
##              icon = "question", type = "okcancel")


###################################################
### code chunk number 118: tkwait
###################################################
msg <- "We care ..."
dialog <- tktoplevel(); tkwm.withdraw(dialog)
tkwm.overrideredirect(dialog, TRUE)   # no decoration
frame <- ttkframe(dialog, padding = 5)
tkpack(frame, expand = TRUE, fill = "both")
tkpack(ttklabel(frame, text = msg), pady = 5)


###################################################
### code chunk number 119: waitVariable
###################################################
flag <- tclVar("")
tkpack(ttkbutton(frame, text="dismiss", command=function() {
  tkgrab.release(dialog)
  tclvalue(flag) <- "Destroy"
}))


###################################################
### code chunk number 120: Dialogs.Rnw:72-74 (eval = FALSE)
###################################################
## tkwm.deiconify(dialog)
## tkwait.variable(flag)


###################################################
### code chunk number 121: getOpen (eval = FALSE)
###################################################
## tkgetOpenFile(filetypes = paste(
##                 "{{jpeg files} {.jpg .jpeg} }",
##                 "{{png files} {.png}}",
##                 "{{All files} {*}}", sep = " ")) # needs space


###################################################
### code chunk number 122: Dialogs.Rnw:121-147
###################################################
window <- tktoplevel()
tkwm.title(window, "File menu example")
menu_bar <- tkmenu(window)
tkconfigure(window, menu = menu_bar)
file_menu <- tkmenu(menu_bar)
tkadd(menu_bar, "cascade", label="File", menu = file_menu)

tkadd(file_menu,"command", label = "Source file...",
      command =  function() {
        file_name <- tkgetOpenFile(filetypes=
                        "{{R files} {.R}} {{All files} *}")
        if(file.exists(file_name <- as.character(file_name)))
           source(tclvalue(file_name))
      })
tkadd(file_menu, "command", label = "Save workspace as...",
      command = function() {
        file_name <- tkgetSaveFile(defaultextension = "Rsave")
        if(nchar(fname <- as.character(file_name)))
          save.image(file = file_name)
      })
tkadd(file_menu, "command", label="Set working directory...",
      command = function() {
        dir_name <- tkchooseDirectory()
        if(nchar(dir_name <- as.character(dir_name)))
          setwd(dir_name)
      })


###################################################
### code chunk number 123: ColorSelection
###################################################
window <- tktoplevel()
tkwm.title(window, "Select a color")
frame <- ttkframe(window, padding = c(3,3,3,12))
tkpack(frame, expand = TRUE, fill = "both")
color_well <- tkcanvas(frame, width = 40, height = 16, 
                      background = "#ee11aa",
                      highlightbackground = "#ababab") 
tkpack(color_well)
tkpack(ttklabel(frame, text = "Click color to change"))
#
tkbind(color_well,"<Button-1>", function(W) {
  color <- tcl("tk_chooseColor", parent = W, 
               title = "Set box color")
  color <- tclvalue(color)
  print(color)
  if(nchar(color))
    tkconfigure(W, background = color)
})


###################################################
### code chunk number 124: Widgets.Rnw:34-37
###################################################
value <- tclVar(TRUE)
tclvalue(value) <- FALSE
tclvalue(value)


###################################################
### code chunk number 125: setup-window
###################################################
window <- tktoplevel(); tkwm.title(window, "Check button example")
frame <- ttkframe(window); tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 126: make-TCL-variables
###################################################
value_var <- tclVar(TRUE)
callback <- function() print(tclvalue(value_var)) # uses global
label_var <- tclVar("check button label")
check_button <- 
  ttkcheckbutton(frame, variable = value_var, 
                 textvariable = label_var, command = callback)
tkpack(check_button)


###################################################
### code chunk number 127: Widgets.Rnw:99-100
###################################################
tkconfigure(check_button, style = "Toolbutton")


###################################################
### code chunk number 128: Widgets.Rnw:112-113
###################################################
tkcget(check_button, "variable" = NULL)


###################################################
### code chunk number 129: ourTclVar
###################################################
our_tcl_var <- function(...) {
  var <- tclVar(...)
  .TkRoot$env[[as.character(var)]] <- var
  var
}
## lookup function
get_tcl_var_by_id <- function(id) {
  .TkRoot$env[[as.character(id)]]
}


###################################################
### code chunk number 130: our_tcl_varExample
###################################################
window <- tktoplevel(); tkwm.title(window, "Check button example")
frame <- ttkframe(window); tkpack(frame, expand = TRUE, fill = "both")
value_var <- our_tcl_var(TRUE)


###################################################
### code chunk number 131: Widgets.Rnw:142-146
###################################################
callback <- function(W) {
  id <- tkcget(W, "variable" = NULL)
  print(get_tcl_var_by_id(id))
}


###################################################
### code chunk number 132: Widgets.Rnw:149-155
###################################################
label_var <- tclVar("check button label")
check_button<- 
  ttkcheckbutton(frame, variable = value_var, 
                 textvariable = label_var)
tkbind(check_button, "<Button-1>", callback)
tkpack(check_button)


###################################################
### code chunk number 133: radio-button-1
###################################################
window <- tktoplevel(); tkwm.title(window, "Radio example")
frame <- ttkframe(window, padding = c(3,3,12,12)); tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 134: radio-button-2
###################################################
values <- c("less", "greater", "two.sided")
var <- tclVar(values[3])                # initial value
callback <- function() print(tclvalue(var))
sapply(values, function(i) {
  radio_button <- ttkradiobutton(frame, variable = var, 
                                 text = i, value = i, 
                                 command = callback)
  tkpack(radio_button, side = "top", anchor = "w")
})


###################################################
### code chunk number 135: entryExample
###################################################
window <- tktoplevel()
tkwm.title(window, "Entry widget test")
frame <- ttkframe(window, padding = c(3,3,12,12)); tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 136: entryExampleDef
###################################################
txt_var <- tclVar("initial value")
entry <- ttkentry(window, textvariable = txt_var)
tkpack(entry)


###################################################
### code chunk number 137: Widgets.Rnw:228-230
###################################################
tclvalue(txt_var)
tclvalue(txt_var) <- "set value"


###################################################
### code chunk number 138: tkget
###################################################
tkget(entry)


###################################################
### code chunk number 139: tkinsert
###################################################
tkinsert(entry, "end", "new text")


###################################################
### code chunk number 140: Widgets.Rnw:254-255
###################################################
tkdelete(entry, 0, 4)


###################################################
### code chunk number 141: Widgets.Rnw:263-264
###################################################
tkicursor(entry, 0)                         # move to beginning


###################################################
### code chunk number 142: Widgets.Rnw:272-273
###################################################
tkselection.range(entry, 0, "end")


###################################################
### code chunk number 143: InitialMsg
###################################################
## "R5" class for ttk entry with initial message
library(tcltk)


###################################################
### code chunk number 144: ex-tcltk-initial-message.Rnw:16-23
###################################################
setOldClass(c("tkwin", "tclVar"))
TtkEntry <- setRefClass("TtkEntry",
                        fields = list(
                          entry = "tkwin",     # entry
                          tcl_var  = "tclVar", # text variable
                          init_msg = "character"
                          ))


###################################################
### code chunk number 145: init-msg-style
###################################################
.Tcl("ttk::style configure Gray.TEntry -foreground gray") 


###################################################
### code chunk number 146: init-msg-methods
###################################################
TtkEntry$methods(
                 is_init_msg = function() {
                   "Is the init text showing?"
                   as.character(tclvalue(tcl_var)) == init_msg
                 },
                 hide_init_msg = function() {
                   "Hide the initial text"
                   if(is_init_msg()) {
                     tkconfigure(entry, style = "TEntry")
                     set_text("", hide = FALSE)
                   }
                 },
                 show_init_msg = function() {
                   "Show the initial text"
                   tkconfigure(entry, style = "Gray.TEntry")
                   set_text(init_msg, hide = FALSE)
                 })


###################################################
### code chunk number 147: get-set-text
###################################################
TtkEntry$methods(
                 set_text = function(text, hide = TRUE) {
                   "Set text into widget"
                   if(hide) hide_init_msg()
                   tcl_var_local <- tcl_var   # avoid warning
                   tclvalue(tcl_var_local) <- text
                 },
                 get_text = function() {
                   "Get the text value"
                   if(!is_init_msg())
                     as.character(tclvalue(tcl_var))
                   else
                     ""
                 })


###################################################
### code chunk number 148: add-bindings
###################################################
TtkEntry$methods(
                 add_bindings = function() {
                   "Add focus bindings to make this work"
                   tkbind(entry, "<FocusIn>", hide_init_msg)
                   tkbind(entry, "<FocusOut>", function() {
                     if(nchar(get_text()) == 0)
                       show_init_msg()
                   })
                 })


###################################################
### code chunk number 149: ex-tcltk-initial-message.Rnw:96-109
###################################################
TtkEntry$methods(
   initialize = function(parent, text, init_msg = "", ...) {
     tcl_var <<- tclVar()
     entry <<- ttkentry(parent, textvariable = tcl_var)
     init_msg <<- init_msg
     ##
     if(missing(text))
       show_init_msg()
     else
       set_text(text)
     add_bindings()
     callSuper(...)
   })


###################################################
### code chunk number 150: ex-tcltk-initial-message.Rnw:116-126
###################################################
window <- tktoplevel()
widget <- TtkEntry$new(parent = window, 
                       init_msg = "type value here")
tkpack(widget$entry)
#
button <- ttkbutton(window, text = "focus out onto this", 
               command = function() {
                 print(widget$get_text())
               })
tkpack(button)


###################################################
### code chunk number 151: ex-tcltk-validation.Rnw:6-40
###################################################
## Example of using validation to adjust the date
## in case a user doesn't use desired format

## Docs on validation
## VALIDATION
## The -validate, -validatecommand, and -invalidcommand options are used to enable entry widget validation.
## VALIDATION MODES
## There are two main validation modes: prevalidation, in which the -validatecommand is evaluated prior to each edit and the return value is used to determine whether to accept or reject the change; and revalidation, in which the -validatecommand is evaluated to determine whether the current value is valid.

## The -validate option determines when validation occurs; it may be set to any of the following values:

## none
##     Default. This means validation will only occur when specifically requested by the validate widget command.

## key
##     The entry will be prevalidated prior to each edit (specifically, whenever the insert or delete widget commands are called). If prevalidation fails, the edit is rejected.

## focus
##     The entry is revalidated when the entry receives or loses focus.

## focusin
##     The entry is revalidated when the entry receives focus.

## focusout
##     The entry is revalidated when the entry loses focus.

## all
##     Validation is performed for all above conditions.

## The -invalidcommand is evaluated whenever the -validatecommand returns a false value.

## The -validatecommand and -invalidcommand may modify the entry widget's value via the widget insert or delete commands, or by setting the linked -textvariable. If either does so during prevalidation, then the edit is rejected regardless of the value returned by the -validatecommand.

## If -validatecommand is empty (the default), validation always succeeds.


###################################################
### code chunk number 152: ex-tcltk-validation.Rnw:44-49
###################################################
## test of validation command
## tricky bit is that validation commands must return TRUE or FALSE
## we can do this with tcl("eval","FALSE") or .Tcl("expr false")

library(tcltk)


###################################################
### code chunk number 153: ex-tcltk-validation.Rnw:85-92
###################################################
date_patterns <- c()
for(i in list(c("%m","%d","%Y"),        # U.S. style
              c("%m","%d","%y"))) {
  for(j in c("/","-"," ") )
    date_patterns[length(date_patterns)+1] <- 
      paste(i,sep = "", collapse = j)
}


###################################################
### code chunk number 154: setValidDateCallback
###################################################
is_valid_date <- function(W, P) { # P is the current value
  for(i in date_patterns) {
    date <- try( as.Date(P, format = i), silent = TRUE)
    if(!inherits(date, "try-error") && !is.na(date)) {
      tkconfigure(W, foreground = "black")  # or use style
      tkdelete(W, 0,"end")
      tkinsert(W, 0, format(date, format = "%m/%d/%y"))
      return(tcl("expr","TRUE"))        
    } 
  }
  return(tcl("expr","FALSE"))
}


###################################################
### code chunk number 155: setInvalidCallback
###################################################
indicate_invalid_date <- function(W) {
  tkconfigure(W, foreground = "red")
  tcl("expr", "TRUE")
}


###################################################
### code chunk number 156: notShown
###################################################
## A simple GUI to show the entry widget.
window <- tktoplevel(); tkwm.title(window, "Validation of date example")
frame <- ttkframe(window, padding = c(3,3,12,12)); tkpack(frame, expand = TRUE, fill = "both")
tkpack(ttklabel(frame, text = "Enter a date (mm/dd/yy):"), side = "left", padx = 2)


###################################################
### code chunk number 157: entryWidgetWithValidation
###################################################
entry <- ttkentry(frame, validate = "focusout",
                  validatecommand = is_valid_date,
                  invalidcommand = indicate_invalid_date)
button <- ttkbutton(frame, text = "click")  # focus target
sapply(list(entry, button), tkpack, side = "left", padx = 2)


###################################################
### code chunk number 158: Widgets.Rnw:327-330
###################################################
window <- tktoplevel(); tkwm.title(window, "Combo box example")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 159: Widgets.Rnw:333-335
###################################################
values <- state.name
var <- tclVar(values[1])              # initial value


###################################################
### code chunk number 160: Widgets.Rnw:339-345
###################################################
combo_box <- ttkcombobox(frame,
                         values = values,
                         textvariable = var,
                         state = "normal",     # or "readonly"
                         justify = "left")
tkpack(combo_box)


###################################################
### code chunk number 161: Combobox-set-values
###################################################
tkconfigure(combo_box, values = tolower(values))


###################################################
### code chunk number 162: combobox-set-length-1 (eval = FALSE)
###################################################
## tkconfigure(combo_box, values = as.tclObj("New York"))


###################################################
### code chunk number 163: Combobox-set
###################################################
tclvalue(var) <- values[2]              # using tcl variable
tkset(combo_box, values[4])             # by value
tcl(combo_box, "current", 4)            # by index


###################################################
### code chunk number 164: Combobox-get
###################################################
tclvalue(var)                           # TCL variable
tkget(combo_box)                        # get subcommand
tcl(combo_box, "current")               # 0-based index


###################################################
### code chunk number 165: combobox-selection-binding
###################################################
tkbind(combo_box, "<<ComboboxSelected>>", function() {
  print(tclvalue(var))
})


###################################################
### code chunk number 166: Combobox-binding-to-return
###################################################
tkbind(combo_box, "<Return>", function(W) {
  val <- tkget(W)
  cat(as.character(val), "\n")
})


###################################################
### code chunk number 167: ttkscale
###################################################
ttkscale <- function(parent, ...) 
  tkwidget(parent, "ttk::scale", ...)


###################################################
### code chunk number 168: ttksliderclass
###################################################
Slider <-
  setRefClass("TtkSlider",
     fields = c("frame", "widget", "var", "x", "FUN"),
     methods = list(
       initialize = function(parent, x, ...) {
         initFields(x = x, var = tclVar(1),
                    FUN = NULL, frame = ttkframe(parent))
         widget <<- ttkscale(frame, from = 1, to = length(x),
                       variable = var, orient = "horizontal")
         ## For this widget, the callback is passed a value 
         ## which we ignore here
         tkconfigure(widget, command = function(...) {
           if(is.function(FUN)) FUN(.self)
         })
         layout_gui()
         callSuper(...)
       },
       layout_gui = function() {         
         tkgrid(widget, row = 0, column = 0, columnspan = 3, 
                sticky = "we")
         tkgrid(ttklabel(frame, text = x[1]), 
                row = 1, column = 0)
         tkgrid(ttklabel(frame, text = x[length(x)]), 
                row = 1, column = 2)
         tkgrid.columnconfigure(frame, 1, weight = 1)
       },
       add_callback = function(FUN) FUN <<- FUN,
       get_value = function() x[as.numeric(tclvalue(var))],
       set_value = function(value) {
         "Set value. Value must be in x"
         ind <- match(value, x)
         if(!is.na(ind)) {
           var_local <- var
           tclvalue(var_local) <- ind
         }
       }
       ))


###################################################
### code chunk number 169: Widgets.Rnw:496-506
###################################################
window <- tktoplevel()
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")
x <- seq(0,1,by = 0.05)
##
slider <- Slider$new(parent = window, x = x)
tkpack(slider$frame, expand = TRUE, fill = "x", anchor = "n")
##
slider$set_value(0.5)
print(slider$get_value())


###################################################
### code chunk number 170: use-slider-command
###################################################
slider$add_callback(function(obj) print(obj$get_value()))


###################################################
### code chunk number 171: tkspinbox
###################################################
tkspinbox <- function(parent, ...) 
    tkwidget(parent, "tk::spinbox", ...)


###################################################
### code chunk number 172: Widgets.Rnw:553-555
###################################################
window <- tktoplevel()
spin_box <- tkspinbox(window, values = state.name, wrap=TRUE)


###################################################
### code chunk number 173: Widgets.Rnw:559-560
###################################################
spin_box1 <- tkspinbox(window, from=1, to = 10, increment = 1)


###################################################
### code chunk number 174: Widgets.Rnw:563-565
###################################################
tkpack(spin_box)
tkpack(spin_box1)


###################################################
### code chunk number 175: ex-tcltk-t-test.Rnw:16-19
###################################################
## t.test dialog
## using basic widgets -- no entry widgets yet
library(tcltk)


###################################################
### code chunk number 176: ex-tcltk-t-test.Rnw:22-40
###################################################
## helper functions
## not shown

get_numeric_vars <- function(DF) {
  if(missing(DF))
    return(c(""))
  ProgGUIinR:::find_vars(DF, is.numeric)
}
get_two_level_factor <- function(DF) {
  if(missing(DF))
    return(c(""))
  nms <- names(DF)
  ind <- sapply(DF, function(i) length(levels(as.factor(i))) == 2)
  if(length(ind) > 0)
    nms[ind]
  else
    c("")
}


###################################################
### code chunk number 177: dataModel
###################################################
e <- new.env()
e$x <- tclVar(""); e$f <- tclVar(""); e$data <- tclVar("")
e$mu <- tclVar(0); e$alternative <- tclVar("two.sided")
e$conf.level <- tclVar(95); e$var.equal <- tclVar(FALSE)


###################################################
### code chunk number 178: ex-tcltk-t-test.Rnw:58-88
###################################################
## We don't show the function runTTest.
## It is a bit long, as care must be taken as it isn't clear if a formula should be used.  
runTTest <- function() {
  l <- lapply(e, tclvalue)
  
  ## ugly function to run t test
  if(l$data == "" || l$x == "")
    return()

  l$data <- get(l$data, envir = .GlobalEnv)

  if(l$f != "") {
    l$formula <- formula(paste(l$x,l$f, sep = "~"))
    l$x <- l$f <- NULL
    l$mu <- NULL
    l$var.equal <- as.logical(as.numeric(l$var.equal))

    TTest <- stats:::t.test.formula
  } else {
    l$x <- l$data[, l$x]
    l$f <- NULL
    l$mu = as.numeric(l$mu)
    l$var.equal <- NULL
    
    TTest <- stats:::t.test.default
  }
  l$conf.level <- as.numeric(l$conf.level)/100
  ret <- capture.output(do.call("TTest", l))
  cat(paste(ret, collapse = "\n"))
}


###################################################
### code chunk number 179: notShown
###################################################
### GUI Our standard setup
window <- tktoplevel()
tkwm.title(window, "t-test Dialog")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 180: layout
###################################################
label_frame <- ttklabelframe(frame, text = "t.test()", 
                             padding = 10)
tkpack(label_frame, expand = TRUE, fill = "both", 
       padx = 5, pady = 5)


###################################################
### code chunk number 181: ex-tcltk-t-test.Rnw:115-119
###################################################
tkgrid.columnconfigure(label_frame, 0, weight = 1)
tkgrid.columnconfigure(label_frame, 1, weight = 10)
tkgrid.columnconfigure(label_frame, 2, weight = 1)
tkgrid.columnconfigure(label_frame, 1, weight = 10)


###################################################
### code chunk number 182: tkgridHelper
###################################################
put_label <- function(parent, text, row, column) {
  label <- ttklabel(parent, text = text)
  tkgrid(label, row = row, column = column, sticky = "e")
}


###################################################
### code chunk number 183: readonly
###################################################
put_label(label_frame, "data:",0,0)
data_combo <- ttkcombobox(label_frame, state = "readonly", 
                         values = ProgGUIinR:::avail_dfs(), 
                         textvariable = e$data)
tkgrid(data_combo, row = 0, column = 1, sticky="ew", padx = 2)
tkfocus(data_combo)                      # give focus


###################################################
### code chunk number 184: notShown
###################################################
## not shown
put_label(label_frame, "x:",1,0)
x_combo <-  ttkcombobox(label_frame, 
                       values = get_numeric_vars(), textvariable = e$x)
tkgrid(x_combo, row = 1, column = 1, sticky = "ew", padx = 2)


###################################################
### code chunk number 185: notShown
###################################################
## not shown
put_label(label_frame, "~ f:",1,2)
factor_combo <-  ttkcombobox(label_frame, values = get_two_level_factor(), textvariable = e$f)
tkgrid(factor_combo, row = 1, column = 3, sticky = "ew", padx = 2)


###################################################
### code chunk number 186: mu
###################################################
put_label(label_frame, "mu:", 2, 0)
mu_combo <-  ttkentry(label_frame,  textvariable = e$mu)
tkgrid(mu_combo, row = 2, column = 1, sticky = "ew", padx = 2)


###################################################
### code chunk number 187: ex-tcltk-t-test.Rnw:172-174
###################################################
ttkscale <- function(parent, ...) tkwidget(parent, "ttk::scale", ...)
tkspinbox <- function(parent, ...) tkwidget(parent, "tk::spinbox", ...)


###################################################
### code chunk number 188: ex-tcltk-t-test.Rnw:176-185
###################################################
put_label(label_frame, "alternative:", 3, 0)
rb_frame <- ttkframe(label_frame)
sapply(c("two.sided","less","greater"), function(i) {
  radio_button <- 
    ttkradiobutton(rb_frame, variable = e$alternative, 
                   text = i, value = i)
  tkpack(radio_button, side = "left")
})
tkgrid(rb_frame, row = 3, column = 1, sticky = "ew", padx = 2)


###################################################
### code chunk number 189: ex-tcltk-t-test.Rnw:191-206
###################################################
put_label(label_frame, "conf.level:", 3, 2)
conf_level_frame <- ttkframe(label_frame)
tkgrid(conf_level_frame, row = 3, column = 3, columnspan = 2, 
       sticky = "ew", padx = 2)
##
conf_level_scale <- ttkscale(conf_level_frame, 
                     from = 75, to = 100,  
                     variable = e$conf.level)
conf_level_spin <- tkspinbox(conf_level_frame, 
                     from = 75, to = 100, increment = 1, 
                     textvariable = e$conf.level, width = 5)
##
tkpack(conf_level_scale, expand = TRUE, fill = "y",
       side = "left")
tkpack(conf_level_spin, side = "left")


###################################################
### code chunk number 190: ex-tcltk-t-test.Rnw:210-215
###################################################
put_label(label_frame, "var.equal:", 4, 0)
var_equal_check <- 
  ttkcheckbutton(label_frame, variable = e$var.equal)
tkgrid(var_equal_check, row = 4, column = 1, stick = "w", 
       padx = 2)


###################################################
### code chunk number 191: ok-cancel
###################################################
button_frame <- ttkframe(frame)
cancel_button <- ttkbutton(button_frame, text = "cancel")
ok_button <- ttkbutton(button_frame, text = "ok")
#
tkpack(button_frame, fill = "x", padx = 5, pady = 5)
tkpack(ttklabel(button_frame, text = " "), expand = TRUE,
       fill = "y", side = "left")               # add a spring
sapply(list(cancel_button, ok_button), tkpack, 
       side = "left", padx = 6)


###################################################
### code chunk number 192: ex-tcltk-t-test.Rnw:237-241
###################################################
tkconfigure(ok_button, command = runTTest)
tkconfigure(cancel_button, 
            command = function() tkdestroy(window))
tkbind("TButton", "<Return>", function(W) tcl(W, "invoke"))


###################################################
### code chunk number 193: ex-tcltk-t-test.Rnw:251-297
###################################################
update_ui <- function() {
  dfName <- tclvalue(e$data)
  curDfs <- ProgGUIinR:::avail_dfs()
  tkconfigure(data_combo, values = curDfs)
  if(!dfName %in% curDfs) {
    dfName <- ""
    tclvalue(e$data) <- ""
  }

  if(dfName == "") {
    ## 3 ways to disable needed!!
    x <- list(x_combo, factor_combo, mu_combo,  
              conf_level_scale, var_equal_check, ok_button)
    sapply(x, function(W) tcl(W, "state", "disabled"))
    sapply(as.character(tkwinfo("children", rb_frame)), 
           function(W) tcl(W, "state", "disabled"))
    tkconfigure(conf_level_spin, state = "disabled")
  } else {
    ## enable univariate, ok
    sapply(list(x_combo,mu_combo,conf_level_scale,ok_button),
           function(W) tcl(W, "state", "!disabled"))
    sapply(as.character(tkwinfo("children", rb_frame)), 
           function(W) tcl(W, "state", "!disabled"))
    tkconfigure(conf_level_spin, state = "normal")
    
    DF <- get(dfName, envir = .GlobalEnv)
    numVars <- get_numeric_vars(DF)
    tkconfigure(x_combo, values = numVars)
    if(! tclvalue(e$x) %in% numVars)
      tclvalue(e$x) <- ""

    ## bivariate
    avail_factors <- get_two_level_factor(DF)
    sapply(list(factor_combo, var_equal_check),
           function(W) {
             val <- if(length(avail_factors)) "!" else ""
             tcl(W, "state", sprintf("%sdisabled", val))
           })
    tkconfigure(factor_combo, values = avail_factors)
    if(!tclvalue(e$f) %in% avail_factors)
      tclvalue(e$f) <- ""
      
         }
}
update_ui()
tkbind(data_combo, "<<ComboboxSelected>>", update_ui)


###################################################
### code chunk number 194: digest
###################################################
require(digest)
create_function <- function() {
  .DFs <- digest(ProgGUIinR:::avail_dfs())
  f <- function(...) {
    if((val <- digest(ProgGUIinR:::avail_dfs())) != .DFs) {
      .DFs <<- val
      update_ui()
    }
    return(TRUE)
  }
}


###################################################
### code chunk number 195: taskcallback (eval = FALSE)
###################################################
## id <- addTaskCallback(create_function())


###################################################
### code chunk number 196: scrollbar-example
###################################################
library(tcltk)
window <- tktoplevel()
tkwm.title(window, "Scroll bar example")
parent <- ttkframe(window)
tkpack(parent, expand = TRUE, fill = "both")
widget <- tktext(parent)


###################################################
### code chunk number 197: ScrollableWidgets.Rnw:26-30
###################################################
xscr <- ttkscrollbar(parent, orient = "horizontal",
                 command = function(...) tkxview(widget, ...))
yscr <- ttkscrollbar(parent, orient = "vertical",
                 command = function(...) tkyview(widget, ...))


###################################################
### code chunk number 198: ScrollableWidgets.Rnw:40-43
###################################################
tkconfigure(widget,
            xscrollcommand = function(...) tkset(xscr,...),
            yscrollcommand = function(...) tkset(yscr,...))


###################################################
### code chunk number 199: ScrollableWidgets.Rnw:50-55
###################################################
tkgrid(widget, row = 0, column = 0, sticky = "news")
tkgrid(yscr,   row = 0, column = 1, sticky = "ns")
tkgrid(xscr,   row = 1, column = 0, sticky = "ew")
tkgrid.columnconfigure(parent, 0, weight = 1)
tkgrid.rowconfigure(parent, 0, weight = 1)


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


###################################################
### code chunk number 213: makeGUI
###################################################
w <- tktoplevel(); tkwm.title(w, "Text buffer example")
f <- ttkframe(w, padding = c(3,3,12,12))
tkpack(f, expand = TRUE, fill = "both")
txt <- tktext(f, width = 80, height = 24)   # default size
addScrollbars(f, txt)


###################################################
### code chunk number 214: ex-tcltk-text.Rnw:24-29
###################################################
tktag.configure(txt, "commandTag", foreground = "blue", 
                font = "courier 12 italic")
tktag.configure(txt, "outputTag", font = "courier 12")
tktag.configure(txt, "errorTag", foreground = "red", 
                font = "courier 12 bold")


###################################################
### code chunk number 215: ex-tcltk-text.Rnw:36-60
###################################################
eval_cmd_chunk <- function(txt, cmds) {
  
  cmd_chunks <- try(parse(text = cmds), silent = TRUE)
  if(inherits(cmd_chunks,"try-error")) {
    tkinsert(t, "end", "Error", "errorTag") # add markup tag
  }

  for(cmd in cmd_chunks) {
    cutoff <- 0.75 * getOption("width")
    dcmd <- deparse(cmd, width.cutoff = cutoff)
    command <- 
      paste(getOption("prompt"),
            paste(dcmd, collapse = paste("\n", 
                          getOption("continue"), sep = "")),
            sep = "", collapse = "")
    tkinsert(txt, "end", command, "commandTag")
    tkinsert(txt, "end","\n")
    ## output, should check for errors in eval!
    output <- capture.output(eval(cmd, envir = .GlobalEnv))
    output <- paste(output, collapse = "\n")
    tkinsert(txt, "end", output, "outputTag")
    tkinsert(txt, "end","\n")
  }
}


###################################################
### code chunk number 216: ex-tcltk-text.Rnw:63-80
###################################################
## function to add scrollbars to a widget
addScrollbars <- function(parent, widget) {
  xscr <- ttkscrollbar(parent, orient = "horizontal",
                       command = function(...) tkxview(widget, ...))
  yscr <- ttkscrollbar(parent, orient = "vertical",
                       command = function(...) tkyview(widget, ...))

  tkconfigure(widget,
              xscrollcommand = function(...) tkset(xscr,...),
              yscrollcommand = function(...) tkset(yscr,...))

  tkgrid(widget, row = 0, column = 0, sticky = "news")
  tkgrid(yscr,row = 0,column = 1, sticky = "ns")
  tkgrid(xscr, row = 1, column = 0, sticky = "ew")
  tkgrid.columnconfigure(parent, 0, weight = 1)
  tkgrid.rowconfigure(parent, 0, weight = 1)
}


###################################################
### code chunk number 217: TryIt
###################################################
eval_cmd_chunk(txt, "2 + 2; lm(mpg ~ wt, data = mtcars)")


###################################################
### code chunk number 218: tearoff
###################################################
tcl("option","add","*tearOff", 0)    # disable tearoff menus


###################################################
### code chunk number 219: testIfMac
###################################################
using_Mac <- function()  
  as.character(tcl("tk", "windowingsystem")) == "aqua"


###################################################
### code chunk number 220: ScrollableWidgets.Rnw:371-377
###################################################
window <- tktoplevel()
tkwm.title(window, "Pop-Up menu example")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")

button <- ttkbutton(frame, text = "Click me for pop-up")


###################################################
### code chunk number 221: ScrollableWidgets.Rnw:380-387
###################################################
doPopup <- function(X, Y) tkpopup(pmb, X, Y) # define callback
if (using_Mac()) {
  tkbind(button, "<Button-2>", doPopup)      # right click
  tkbind(button, "<Control-1>", doPopup)     # Control + click
} else {
  tkbind(button, "<Button-3>", doPopup)
}


###################################################
### code chunk number 222: ex-tcltk-menu.Rnw:1-34
###################################################
## Not shown
library(tcltk)

## Helper functions
using_Mac <- function()  
  as.character(tcl("tk", "windowingsystem")) == "aqua"

addScrollbars <- function(parent, widget,type = c("both", "x", "y")) {
  if(any(type %in% c("both","x"))) {
    xscr <- ttkscrollbar(parent, orient = "horizontal",
                         command = function(...) tkxview(widget, ...))
    tkconfigure(widget,
                xscrollcommand = function(...) tkset(xscr,...))
  }

  if(any(type %in% c("both","y"))) {
     yscr <- ttkscrollbar(parent, orient = "vertical",
                          command = function(...) tkyview(widget, ...))
     tkconfigure(widget,
                 yscrollcommand = function(...) tkset(yscr,...))
   }

  ## place in grid     
  tkgrid(widget, row = 0, column = 0, sticky = "news")
  if(any(type %in% c("both", "x"))) {
    tkgrid(xscr, row = 1, column = 0, sticky = "ew")
    tkgrid.columnconfigure(parent, 0, weight = 1)
  }
  if(any(type %in% c("both", "y"))) {
    tkgrid(yscr,row = 0,column = 1, sticky = "ns")
    tkgrid.rowconfigure(parent, 0, weight = 1)
  }
}


###################################################
### code chunk number 223: ex-tcltk-menu.Rnw:40-44
###################################################
library(svMisc)                         # for some helpers
showCmd <- function(cmd) {
  writeLines(captureAll(parseText(cmd)))
}


###################################################
### code chunk number 224: ex-tcltk-menu.Rnw:49-55
###################################################
window <- tktoplevel()
tkwm.title(window, "Simple code editor")
frame <- ttkframe(window, padding = c(3,3,12,12)) 
tkpack(frame, expand = TRUE, fill = "both")
text_buffer <- tktext(frame, undo = TRUE)
addScrollbars(frame, text_buffer)


###################################################
### code chunk number 225: ex-tcltk-menu.Rnw:61-69
###################################################
menu_bar <- tkmenu(window)
tkconfigure(window, menu = menu_bar)
#
file_menu <- tkmenu(menu_bar)
tkadd(menu_bar, "cascade", label = "File", menu = file_menu)
#
edit_menu <- tkmenu(menu_bar)
tkadd(menu_bar, "cascade", label = "Edit", menu = edit_menu)


###################################################
### code chunk number 226: ex-tcltk-menu.Rnw:75-80
###################################################
tkadd(file_menu, "command", label = "Evaluate buffer",
      command = function() {
        cur_val <- tclvalue(tkget(text_buffer, "1.0", "end"))
        showCmd(cur_val)
      })


###################################################
### code chunk number 227: ex-tcltk-menu.Rnw:84-91
###################################################
tkadd(file_menu, "command", label = "Evaluate selection",
      state = "disabled",
      command =  function() {
        cur_sel <- tclvalue(tkget(text_buffer,
                                  "sel.first", "sel.last"))
        showCmd(cur_sel)
      })


###################################################
### code chunk number 228: addQuit
###################################################
tkadd(file_menu, "separator")
tkadd(file_menu, "command", label = "Quit", 
      command = function() tkdestroy(window))


###################################################
### code chunk number 229: ex-tcltk-menu.Rnw:102-110
###################################################
img <- system.file("images","up.gif", package = "gWidgets")
tkimage.create("photo", "::img::undo", file = img)
tkadd(edit_menu, "command", label = "Undo",
      image = "::img::undo", compound = "left", 
      state = "disabled",
      command = function() tcl(text_buffer, "edit", "undo"))
tkadd(edit_menu, "command", label="Redo", state = "disabled",
      command = function() tcl(text_buffer, "edit", "redo"))


###################################################
### code chunk number 230: ex-tcltk-menu.Rnw:116-125
###################################################
tkbind(text_buffer, "<<Selection>>", function(W) {
  hasSelection <- function(W) {
    ranges <- tclvalue(tcl(W, "tag", "ranges", "sel"))
    length(ranges) > 1 || ranges != ""
  }
  ## configure using an index
  sel_state <- ifelse(hasSelection(W), "normal", "disabled")
  tkentryconfigure(file_menu, 2, state = sel_state)
})


###################################################
### code chunk number 231: ex-tcltk-menu.Rnw:129-136
###################################################
tkbind(text_buffer, "<<Modified>>", function(W) {
  ## not really can_undo/can_redo but nothing suitable
  can_undo <- as.logical(tcl(W,"edit", "modified"))
  undo_state <- ifelse(can_undo, "normal", "disabled")
  sapply(c("Undo", "Redo"), function(i)        # match pattern
         tkentryconfigure(edit_menu, i, state = undo_state)) 
})


###################################################
### code chunk number 232: addKeyboardShortcut
###################################################
if(using_Mac()) {
  tkentryconfigure(edit_menu, "Undo", accelerator="Cmd-z")
  tkbind(window, "<Option-z>", function() {
    tcl(text_buffer, "edit", "undo")
  })
} else {
  tkentryconfigure(edit_menu, "Undo", accelerator="Control-u")
  tkbind(window, "<Control-u>", function() {
    tcl(text_buffer, "edit", "undo")
  })
}


###################################################
### code chunk number 233: definePopup
###################################################
do_popup <- function(W, X, Y) {
  cur <- tclvalue(tkget(W, "current  wordstart", 
                           "current wordend"))
  tcl(W, "tag", "add", "popup", "current  wordstart", 
                                "current wordend")
  possible_vals <- head(completion(cur)[,1, drop=TRUE], n=20)
  if(length(possible_vals) > 1) {
    popup <- tkmenu(text_buffer)       # create menu for popup
    sapply(possible_vals, function(i) {         
      tkadd(popup, "command", label=i, command = function() {
        tcl(W,"replace", "popup.first", "popup.last", i)
      })
    })
    tkpopup(popup, X, Y)
 }}


###################################################
### code chunk number 234: addPopup
###################################################
if (!using_Mac()) {
  tkbind(text_buffer, "<Button-3>", do_popup)
} else {
  tkbind(text_buffer, "<Button-2>", function(W,X,Y) {
    ## UNIX legacy re mouse-2 click for selection copy
    tcl("clipboard","clear",displayof = W) 
    do_popup(W,X,Y)
    })      # right click
  tkbind(text_buffer, "<Control-1>", do_popup) # Ctrl+click
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


###################################################
### code chunk number 242: showChildrenIndex
###################################################
children <- tcl(treeview, "children", "")
(children <- head(as.character(children)))     # as.character
sapply(children, function(i) tclvalue(tkindex(treeview, i)))


###################################################
### code chunk number 243: getValue
###################################################
x <- tcl(treeview, "item", children[1], "-values") # no tkitem
as.character(x)


###################################################
### code chunk number 244: ScrollableWidgets.Rnw:722-723
###################################################
tkselect(treeview, "set", children)


###################################################
### code chunk number 245: ScrollableWidgets.Rnw:727-728
###################################################
tkselect(treeview, "toggle", tcl(treeview, "children", ""))


###################################################
### code chunk number 246: ScrollableWidgets.Rnw:732-733
###################################################
IDs <- as.character(tkselect(treeview))


###################################################
### code chunk number 247: ScrollableWidgets.Rnw:747-752
###################################################
callback_example <- function(W, x, y) {
  col <- as.character(tkidentify(W, "column", x, y))
  row <- as.character(tkidentify(W, "row", x, y))
  ## now do something ...
}


###################################################
### code chunk number 248: ex-tcltk-table
###################################################
library(tcltk)

## helpers
quoteIt <- function(string) {           
  doQuote <- function(x) {
    xx <- strsplit(x, '"', fixe = TRUE)[[1]]
    paste(paste('"', xx, '"', sep = ""), collapse = '\'"\'')
  }
  if(!length(string)) return("")
  has_double_quote <- grep('"',string)
  if(!length(has_double_quote))
    return(paste('"',string,'"',sep = ""))
  if (!length(grep("([$`])", string))) {
    paste("\"", gsub("([\"!\\])", "\\\\\\1", string), 
          "\"", sep = "")
  } else sapply(string, doQuote)
}


## covert a data frame into a character based on
.toCharacter <- function(x,width,...) UseMethod(".toCharacter")
.toCharacter.default <- function(x,width,...) as.character(x)
.toCharacter.integer <- function(x,width,...) {
 if(missing(width)) width <- max(nchar(as.character(x))) + 2  
  format(x, justify = "right", width = width)
}
.toCharacter.numeric <- function(x,width,...) {
  if(missing(width)) width <- max(nchar(as.character(x))) + 2
  format(x,trim = FALSE, width = width, justify = "right")
}
.toCharacter.factor <- function(x,width,...) {
  if(missing(width)) width <- max(nchar(as.character(x))) + 2
  .toCharacter(as.character(x),width,...)
}
.toCharacter.logical <- function(x,width,...) {
  if(missing(width)) width <- 7
  format(as.character(x), justify = "centre", width = width)
}
.toCharacter.data.frame <- function(x,width =  10, ...) {
  nms <- dimnames(x)
  DF <- as.data.frame(lapply(x,function(i) .toCharacter(i, width = width)),
                      stringsAsFactors = FALSE)
  dimnames(DF) <- nms
  return(DF)
}

addScrollbars <- function(parent, widget) {
  xscr <- ttkscrollbar(parent, orient = "horizontal",
                       command = function(...) tkxview(widget, ...))
  yscr <- ttkscrollbar(parent, orient = "vertical",
                       command = function(...) tkyview(widget, ...))

  tkconfigure(widget,
              xscrollcommand = function(...) tkset(xscr,...),
              yscrollcommand = function(...) tkset(yscr,...))

  tkgrid(widget, row = 0, column = 0, sticky = "news")
  tkgrid(yscr,row = 0,column = 1, sticky = "ns")
  tkgrid(xscr, row = 1, column = 0, sticky = "ew")
  tkgrid.columnconfigure(parent, 0, weight = 1)
  tkgrid.rowconfigure(parent, 0, weight = 1)
}


###################################################
### code chunk number 249: ex-tcltk-table.Rnw:67-68
###################################################
DF <- getCRANmirrors()[, c(1,2,5,4)]


###################################################
### code chunk number 250: notShown
###################################################
window <- tktoplevel()
tkwm.title(window, "Choose a CRAN mirror")
frame <- ttkframe(window, padding = c(3,3,3,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 251: ex-tcltk-table.Rnw:79-85
###################################################
frame_0 <- ttkframe(frame); tkpack(frame_0, fill = "x")
label <- ttklabel(frame_0, text = "filter:")
tkpack(label, side = "left")
filter_var <- tclVar("")
filter_entry <- ttkentry(frame_0, textvariable = filter_var)
tkpack(filter_entry, side = "left")


###################################################
### code chunk number 252: makeTreeview
###################################################
frame_1 <- ttkframe(frame)
tkpack(frame_1, expand = TRUE, fill = "both")
treeview <- ttktreeview(frame_1, columns = 1:ncol(DF), 
                  displaycolumns = 1:(ncol(DF) - 1), 
                  show = "headings",     # not "tree" 
                  selectmode = "browse") # single selection
addScrollbars(frame_1, treeview)


###################################################
### code chunk number 253: configureColumns
###################################################
widths <- c(100, 75, 400)            # hard coded
nms <- names(DF)
for(i in 1:3) {
  tcl(treeview, "heading", i, text = nms[i])
  tcl(treeview, "column", i, width = widths[i], 
      stretch = TRUE, anchor = "w")
}


###################################################
### code chunk number 254: ex-tcltk-table.Rnw:125-136
###################################################
fillTable <- function(treeview, DF) {
  children <- as.character(tcl(treeview, "children", ""))
  for(i in children) 
    tcl(treeview, "delete", i)                 # out with old
  shade <- c("none", "gray")
  for(i in seq_len(nrow(DF))) 
    tcl(treeview, "insert", "", "end", tag = shade[i %% 2], 
        text = "",  
        values = unlist(DF[i,]))               # in with new
  tktag.configure(treeview, "gray", background = "gray95")
}


###################################################
### code chunk number 255: ex-tcltk-table.Rnw:140-141
###################################################
fillTable(treeview, DF)


###################################################
### code chunk number 256: ex-tcltk-table.Rnw:146-155
###################################################
cur_ind <- 1:nrow(DF)
tkbind(filter_entry, "<KeyRelease>", function(W, K) {
  val <- tclvalue(tkget(W))
  poss_vals <- apply(DF, 1, function(...) 
                    paste(..., collapse = " "))
  ind<- grep(val, poss_vals)
  if(length(ind) == 0) ind <- 1:nrow(DF)
  fillTable(treeview, DF[ind,])
})


###################################################
### code chunk number 257: ex-tcltk-table.Rnw:160-169
###################################################
tkbind(treeview, "<Double-Button-1>", function(W, x, y) {
  sel <- as.character(tcl(W, "identify", "row", x, y))
  vals <- tcl(W, "item", sel, "-values")
  URL <- as.character(vals)[4]          # not tclvalue
  repos <- getOption("repos")
  repos["CRAN"] <- gsub("/$", "", URL[1L])
  options(repos = repos)
  tkwm.withdraw(tkwinfo("toplevel", W))
})


###################################################
### code chunk number 258: ex-tcltk-subset-filter.Rnw:34-35
###################################################
library(tcltk)


###################################################
### code chunk number 259: FilterList
###################################################
setOldClass("tkwin")
setOldClass("tclVar")
FilterList <- setRefClass("FilterList",
                          fields = list(
                            DF = "data.frame",
                            l = "list",
                            id = "ANY",
                            frame = "tkwin" 
                            ))


###################################################
### code chunk number 260: FilterListMethods
###################################################
FilterList$methods(
          setup_gui = function(parent) {
            enc_frame <- ttkframe(parent, padding = 5)
            tkpack(enc_frame, expand = TRUE, fill = "both")
            frame <<- ttkframe(enc_frame)
            button_frame <- ttkframe(enc_frame)
            ## use grid to manage these
            tkgrid(frame, sticky = "news")
            tkgrid(button_frame, sticky = "new")
            tkgrid.rowconfigure(enc_frame, 1, weight = 1)
            tkgrid.columnconfigure(enc_frame, 0, weight = 1)
            ##
            add_button <- 
              ttkbutton(button_frame, text = "Add", 
                        command = function() .self$add())
            preview_button <- 
              ttkbutton(button_frame, text = "Preview", 
                        command = function() .self$preview())
            ##
            sapply(list(add_button, preview_button), tkpack, 
                   side = "left", padx = 5)
          })


###################################################
### code chunk number 261: FilterListInitialize
###################################################
FilterList$methods(
           initialize = function(DF, parent, ...) {
            initFields(DF = DF, l = list(), id = 0L)
            setup_gui(parent)
             callSuper(...)
           })


###################################################
### code chunk number 262: FilterListSelectVariable
###################################################
FilterList$methods(
           select_variable = function() {
             "Return a variable name from the data frame"
             x <- sapply(DF, function(i) class(i)[1])
             m <- cbind(Variables = names(x), Type = x)
             window <- tktoplevel()
             fr <- ttkframe(window, padding = c(3,3,3,12))
             tkpack(fr, expand = TRUE, fill = "both")
             ##
             a <- populate_rectangular_treeview(fr, m)
             tkconfigure(a$frame, width = 300, height = 200)
             tkpack(a$frame, expand = TRUE, fill = "both")
             ## select a value, store in out
             out <- NA
             tkbind(a$tr, "<<TreeviewSelect>>", function(W) {
               sel <- tcl(W, "selection")
               val <- tcl(W, "item", sel, "-values")
               assign("out", as.character(val)[1], 
                      inherits = TRUE)
               tkdestroy(window)
             })
             tkwait.window(window)
             return(out)
           })


###################################################
### code chunk number 263: FilterListAdd
###################################################
FilterList$methods(
           add = function(variable_name, ...) {
             if(missing(variable_name)) 
               variable_name <- select_variable()
             x <- get(variable_name, DF)
             ## new item
             id <<- id + 1
             item <- newFilterItem(x,variable_name, id, .self)
             ## make frame
             enc_frame <- ttkframe(frame)
             tkpack(enc_frame, 
                    expand = TRUE, fill = "both", pady = 2)
             l[[as.character(id)]] <<- list(frame = enc_frame, 
                                            item = item)
             item$make_gui(enc_frame)
           })


###################################################
### code chunk number 264: FilterListRemove
###################################################
FilterList$methods(
           remove=function(id_obj, ...) {
             "Remove. id is character or item object"
             if(!is.character(id_obj))
               id_obj <- id_obj$id
             tkpack.forget(l[[id_obj]]$frame)
             l[[id_obj]] <<- NULL
           })


###################################################
### code chunk number 265: FilterListGetValue
###################################################
FilterList$methods(
           get_value = function() {
             "Return logical value for all filter items"
             if(length(l) == 0)
               return(rep(TRUE, length=nrow(DF)))
             ##
             out <- sapply(l, function(i) i$item$get_value())
             out[is.na(out)] <- FALSE   ## coerce NA to FALSE
             apply(out, 1, all)
           })


###################################################
### code chunk number 266: ex-tcltk-subset-filter.Rnw:227-252
###################################################
FilterList$methods(
           preview = function() {
             "Preview data frame"
             ind <- get_value()
             if(!any(ind)) {
               message("No matches")
               return()
             }
             ## coerce to character
             m <- DF[ind,]
             for(i in seq_along(m)) 
               m[,i] <- as.character(m[,i])
             ##
             window <- tktoplevel()
             fr <- ttkframe(window, padding = c(3,3,3,12))
             tkpack(fr, expand = TRUE, fill = "both")
             a <- populate_rectangular_treeview(fr, m)
             tkconfigure(a$frame, width = 400, height = 300)
             tkpack(a$frame, expand = TRUE, fill = "both")
             ##
             button <- ttkbutton(fr, text = "dismiss", 
                         command=function() tkdestroy(window))
             tkpack(button, anchor = "sw")
             tkwait.window(window)
           })


###################################################
### code chunk number 267: runIt (eval = FALSE)
###################################################
## window <- tktoplevel()
## require(MASS)
## filter_list <- FilterList$new(DF = Cars93, parent = window)


###################################################
### code chunk number 268: newFilterItem
###################################################
newFilterItem <- function(x, nm = deparse(substitute(x)), id, 
                          list_ref) UseMethod("newFilterItem")
newFilterItem.default <- function(x,nm=deparse(substitute(x)), 
                                  id, list_ref) {
  FilterItemNumeric$new(x = x, nm = nm, id = id, 
                        list_ref = list_ref)
}


###################################################
### code chunk number 269: ex-tcltk-subset-filter.Rnw:280-288
###################################################
## not shown
newFilterItem.character <- function(x,  nm = deparse(substitute(x)), id, list_ref) {
  FilterItemCharacter$new(x = x, nm = nm, id = id, list_ref = list_ref)
}

newFilterItem.factor <- function(x,  nm = deparse(substitute(x)), id, list_ref) {
  newFilterItem(as.character(x), nm, id, list_ref)
}


###################################################
### code chunk number 270: FilterItem
###################################################
FilterItem <- setRefClass("FilterItem",
                          fields = list(
                            x = "ANY",
                            nm = "character",
                            id = "character",
                            list_ref = "ANY"
                            ))


###################################################
### code chunk number 271: FilterItemInitialize
###################################################
FilterItem$methods(
           initialize = function(...) {
             initFields(...)
             .self
           },
           get_value = function() {
             "Return logical value of length x"
             stop("Must be subclassed")
           },
           remove = function() list_ref$remove(.self),
           make_gui = function(parent, ...) {
             "Set up GUI, including defining widgets"
             remove_button <- ttkbutton(parent, text="remove",
                                    command = function() {
                                      .self$remove()
                                    })
             tkpack(remove_button, side = "right")
           })


###################################################
### code chunk number 272: FilterItemNumeric
###################################################
FilterItemNumeric <- setRefClass("FilterItemNumeric",
                                 contains = "FilterItem",
                                 fields = list(
                                   ineq_variable = "tclVar",
                                   value_variable = "tclVar"
                                   ))


###################################################
### code chunk number 273: FilterItemNumericGetValue
###################################################
FilterItemNumeric$methods(
      get_value = function() {
        xpr <- paste(nm, tclvalue(ineq_variable), 
                     tclvalue(value_variable))
        eval(parse(text = xpr), 
             envir = list_ref$DF, parent.frame())
      })



###################################################
### code chunk number 274: FilterItemNumericMakeGui
###################################################
FilterItemNumeric$methods(
      make_gui = function(parent) {
        ## standard width for label
        label_width <- max(sapply(names(list_ref$DF), nchar))
        label <- ttklabel(parent, text=nm, width=label_width)
        ## ineq combo
        vals <- c(">=", ">", "==", "!=", "<", "<=")
        ineq_variable <<- tclVar("<=")
        ineq <- ttkcombobox(parent, values = vals, 
                   textvariable = ineq_variable, width = 4)
        ## entry
        value_variable <<- tclVar(max(x, na.rm = TRUE))
        val <- ttkentry(parent, textvariable = value_variable)
        ##
        sapply(list(label, ineq, val), tkpack, side = "left",
               padx = 5)
        callSuper(parent)
      })



###################################################
### code chunk number 275: FilterItemCharacter
###################################################
FilterItemCharacter <- 
  setRefClass("FilterItemCharacter",
              contains = "FilterItem",
              fields = list(
                tr = "tkwin",
                button = "tkwin",
                poss_vals = "character",
                cur_vals = "character"
                ))


###################################################
### code chunk number 276: ex-tcltk-subset-filter.Rnw:443-447
###################################################
FilterItemCharacter$methods(
          get_value = function() {
            x %in% cur_vals
          })


###################################################
### code chunk number 277: sel_by_name
###################################################
sel_by_name <- function(tr, nms) {
  all_ind <- as.character(tcl(tr, "children", ""))
  vals <- sapply(all_ind, function(i) {
    as.character(tcl(tr, "item", i, "-values"))
  })
  ind <- names(vals[vals %in% nms])
  sapply(ind, function(i) tcl(tr, "selection", "add", i))
  sapply(setdiff(all_ind, ind), 
         function(i) tcl(tr, "selection", "remove", i))
}


###################################################
### code chunk number 278: FilterItemShortenPoss
###################################################
FilterItemCharacter$methods(ellipsize = function() {
            tmp <- paste(cur_vals, collapse = ", ")
            if((N <- nchar(tmp)) > 50)
              tmp <- sprintf("%s...%s", substr(tmp, 0, 15),
                             substr(tmp, N-12, N))
            sprintf("%50s", tmp)
          })


###################################################
### code chunk number 279: FilterItemCharacterSelectValuesDialog
###################################################
FilterItemCharacter$methods(
          select_values_dialog = function() {
            window <- tktoplevel()
            fr <- ttkframe(window, padding = c(3,3,12,12))
            tkpack(fr, expand = TRUE, fill = "both")
            tkpack(ttklabel(fr, 
              text = "Select values by extending selection"))
            ## selection
            m <- matrix(poss_vals)
            colnames(m) <- "Values"
            a <- populate_rectangular_treeview(fr, m)
            tkconfigure(a$tr, selectmode = "extended")
            tkconfigure(a$frame, width = 200, height = 300)
            tkpack(a$frame, expand = TRUE, fill = "both")
            
            sel_by_name(a$tr, cur_vals)         # see above
            
            tkbind(a$tr, "<<TreeviewSelect>>", function() {
              ind <- as.character(tcl(a$tr, "selection"))
              cur <- sapply(ind, function(i) {
                as.character(tcl(a$tr, "item", i, "-values"))
              })
              if(length(cur) == 0)
                cur <- character(0)
              cur_vals <<- cur
            })
            ## buttons
            frame_1 <- ttkframe(fr)
            tkpack(frame_1)
            toggle_button <- ttkbutton(frame_1, text="toggle", 
                         command=function() toggle_sel(a$tr))
            set_button <- ttkbutton(frame_1, text = "set", 
                         command=function() tkdestroy(window))
            sapply(list(toggle_button, set_button), tkpack, 
                   side = "left", padx = 5)
            ## make modal
            tkwait.window(window)
            tkconfigure(button, text = ellipsize())
          })


###################################################
### code chunk number 280: FilterItemCharacterMakeGUI
###################################################
FilterItemCharacter$methods(make_gui = function(parent) {
            poss_vals <<- sort(unique(as.character(x)))
            cur_vals <<- poss_vals
            ## label, ineq, val
            l_width <- max(sapply(names(list_ref$DF), nchar))
            label <- ttklabel(parent, text = nm, 
                              width = l_width)
            ##
            in_label <- ttklabel(parent, text = "%in%")
            ##
            button <<- ttkbutton(parent, text = ellipsize(), 
                       command = .self$select_values_dialog)
            ##
            sapply(list(label, in_label), tkpack,
                   side = "left", padx = 5)
            tkpack(button, 
                   expand = TRUE, fill = "x", side = "left")
            callSuper(parent)
          })


###################################################
### code chunk number 281: Helpers
###################################################
## helpers
## We don't show these
addScrollbars <- function(parent, widget) {
  xscr <- ttkscrollbar(parent, orient = "horizontal",
                       command = function(...) tkxview(widget, ...))
  yscr <- ttkscrollbar(parent, orient = "vertical",
                       command = function(...) tkyview(widget, ...))

  tkconfigure(widget,
              xscrollcommand = function(...) tkset(xscr,...),
              yscrollcommand = function(...) tkset(yscr,...))

  tkgrid(widget, row = 0, column = 0, sticky = "news")
  tkgrid(yscr,row = 0,column = 1, sticky = "ns")
  tkgrid(xscr, row = 1, column = 0, sticky = "ew")
  tkgrid.columnconfigure(parent, 0, weight = 1)
  tkgrid.rowconfigure(parent, 0, weight = 1)
}

populate_rectangular_treeview <- function(parent, m) {
  enc_frame <- ttkframe(parent)
  frame <- ttkframe(enc_frame)
  tkpack(frame, expand = TRUE, fill = "both")
  tr <- ttktreeview(frame,
                    columns = seq_len(ncol(m)),
                    show = "headings",
                    selectmode = "browse"
                    )
  addScrollbars(frame, tr)
  tkpack.propagate(enc_frame, FALSE)


  ## headings,widths
  charWidth <- as.integer(tclvalue(tcl("font","measure","TkTextFont","0")))
  sapply(seq_len(ncol(m)), function(i) {
    tcl(tr, "heading", i, text = colnames(m)[i])
    tcl(tr, "column", i, width = 10 + charWidth*max(apply(m, 2, nchar)))
  })
  tcl(tr, "column", ncol(m), stretch = TRUE)
  
  ## values
  if(ncol(m) == 1)  m <- as.matrix(paste("{",m ,"}", sep=""))
  apply(m, 1, function(vals) 
    tcl(tr, "insert", "", "end", values = vals)
        )
  return(list(tr = tr, frame = enc_frame))
}

   

cur_sel <- function(tr) {
  ind <- as.character(tcl(tr, "selection"))
  sapply(ind, function(i) {
    as.character(tcl(tr, "item", i, "-values"))
  })
}

toggle_sel <- function(tr) {
  children <- as.character(tcl(tr, "children", ""))
  tcl(tr, "selection", "toggle", children) 
}


###################################################
### code chunk number 282: ex-tcltk-subset-filter.Rnw:618-622
###################################################
## Call when all is said and done
window <- tktoplevel()
require(MASS)
filter_list <- FilterList$new(DF = Cars93, parent = window)


###################################################
### code chunk number 283: ex-tcltk-tree.Rnw:3-24
###################################################
## not shown
## load in tcltk
library(tcltk)

## helper function to add scrollbars
addScrollbars <- function(parent, widget) {
  xscr <- ttkscrollbar(parent, orient = "horizontal",
                       command = function(...) tkxview(widget, ...))
  yscr <- ttkscrollbar(parent, orient = "vertical",
                       command = function(...) tkyview(widget, ...))

  tkconfigure(widget,
              xscrollcommand = function(...) tkset(xscr,...),
              yscrollcommand = function(...) tkset(yscr,...))

  tkgrid(widget, row = 0, column = 0, sticky = "news")
  tkgrid(yscr,row = 0,column = 1, sticky = "ns")
  tkgrid(xscr, row = 1, column = 0, sticky = "ew")
  tkgrid.columnconfigure(parent, 0, weight = 1)
  tkgrid.rowconfigure(parent, 0, weight = 1)
}


###################################################
### code chunk number 284: ex-tcltk-tree.Rnw:45-50
###################################################
library(XML)
file_name <- "http://www.omegahat.org/RSXML/shortIntro.html"
doc <- htmlTreeParse(file_name, useInternalNodes = TRUE, 
                     error = function(...) {})
root <- xmlRoot(doc)


###################################################
### code chunk number 285: notShown
###################################################
window <- tktoplevel()
tkwm.title(window, "Treeview example with XML")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 286: ex-tcltk-tree.Rnw:61-64
###################################################
treeview <- ttktreeview(frame, displaycolumns = "#all", 
                        columns = 1)
addScrollbars(frame, treeview)                    


###################################################
### code chunk number 287: columnConfiguration
###################################################
tcl(treeview, "heading", "#0", text = "Name")
tcl(treeview, "column",  "#0", minwidth = 20)
tcl(treeview, "heading",  1,   text = "value")
tcl(treeview, "column",   1,   minwidth = 20)


###################################################
### code chunk number 288: ex-tcltk-tree.Rnw:75-95
###################################################
## http://www.omegahat.org/RSXML/shortIntro.html
## xmlChildren gives children
## xmlName gives name of node
## xmlValue gives values stored in node -- text

## issue with quoting values of tree. This taken from shout
quoteIt <- function(string) {           
  doQuote <- function(x) {
    xx <- strsplit(x, '"', fixe = TRUE)[[1]]
    paste(paste('"', xx, '"', sep = ""), collapse = '\'"\'')
  }
  if(!length(string)) return("")
  has_double_quote <- grep('"',string)
  if(!length(has_double_quote))
    return(paste('"',string,'"',sep = ""))
  if (!length(grep("([$`])", string))) {
    paste("\"", gsub("([\"!\\])", "\\\\\\1", string), 
          "\"", sep = "")
  } else sapply(string, doQuote)
}


###################################################
### code chunk number 289: ex-tcltk-tree.Rnw:104-118
###################################################
insertChild <- function(treeview, node, parent = "") {
  l <- list(treeview, "insert", parent, "end", 
            text = xmlName(node))
  children <- xmlChildren(node)
  if(length(children) == 0) {         # add in values
    values <- paste(xmlValue(node), sep = " ", collapse = " ")
    l$values <- as.tclObj(values)     # avoid split on spaces
  }
  tree_path <- do.call("tcl", l)

  if(length(children))                          # recurse
    for(i in children) insertChild(treeview, i, tree_path)
}
insertChild(treeview, root)


###################################################
### code chunk number 290: ex-tcltk-tree.Rnw:131-133
###################################################
.selected_id <- ""                               # globals
.dragging <- FALSE


###################################################
### code chunk number 291: ex-tcltk-tree.Rnw:137-140
###################################################
tkbind(treeview, "<Button-1>", function(W, x, y) {
  .selected_id <<- as.character(tcl(W, "identify","row", x, y))
})  


###################################################
### code chunk number 292: ex-tcltk-tree.Rnw:145-149
###################################################
tkbind(treeview, "<B1-Motion>", function(W, x, y, X, Y) {
  tkconfigure(W, cursor = "diamond_cross")
  .dragging <<-TRUE
})


###################################################
### code chunk number 293: ex-tcltk-tree.Rnw:159-169
###################################################
tkbind(treeview, "<ButtonRelease-1>", function(W, x, y, X, Y) {
  if(.dragging && .selected_id != "") {
    w <- tkwinfo("containing", X, Y)
    if(as.character(w) == as.character(W)) {
      dropID <- as.character(tcl(W, "identify","row", x, y))
      try(tkmove(W, .selected_id, dropID, "0"), silent = TRUE)
    }
  }
  .dragging <<- FALSE; .selected_id <<- "" # reset
})


###################################################
### code chunk number 294: walkTreeReturnAList
###################################################
tree_to_list <- function(treeview) {
  l <- list()
  walk_tree <- function(child, l) {
    l$name <- tclvalue(tcl(treeview,"item", child, "-text"))
    l$value <- as.character(tcl(treeview,"item", child,
                                "-values"))
    children <- as.character(tcl(treeview, "children", child)) 
    if(length(children)) {
      l$children <- list()
      for(i in children) 
        l$children[[i]] <- walk_tree(i, list()) # recurse
    }
    return(l)
  }
  walk_tree("", l)
}



###################################################
### code chunk number 295: ex-tcltk-scrollable-frame.Rnw:1-20
###################################################
## This is also an example of using a canvas to make a scrollable box container
## cf http://mail.python.org/pipermail/python-list/1999-June/005180.html

library(tcltk)
addScrollbars <- function(parent, widget) {
  xscr <- ttkscrollbar(parent, orient = "horizontal",
                       command = function(...) tkxview(widget, ...))
  tkconfigure(widget, xscrollcommand = function(...) tkset(xscr,...))

  yscr <- ttkscrollbar(parent, command = function(...) tkyview(widget,...))
  tkconfigure(widget, yscrollcommand = function(...) tkset(yscr,...))
  
  ## Pack into a grid, from tkFAQ 10.1
  tkgrid(widget,row = 0,column = 0, sticky = "news")
  tkgrid(xscr,row = 1,column = 0, sticky = "ew")
  tkgrid(yscr,row = 0,column = 1, sticky = "ns")
  tkgrid.columnconfigure(parent, 0, weight = 1)
  tkgrid.rowconfigure(parent, 0, weight = 1)
}


###################################################
### code chunk number 296: ex-tcltk-scrollable-frame.Rnw:39-64
###################################################
scrollable_frame <- function(parent, width=300, height=300) {
  canvas_widget <- 
    tkcanvas(parent,
             borderwidth = 0, highlightthickness = 0,
             width = width, height = height)
  addScrollbars(parent, canvas_widget)
  #
  frame <- ttkframe(canvas_widget, padding = c(0,0,0,0))
  frame_id <- tkcreate(canvas_widget, "window", 0, 0, 
                       anchor = "nw", window = frame)
  tkitemconfigure(canvas_widget, frame_id, width = width)
  ## update scroll region
  tkbind(frame, "<Configure>", function() {  
    bbox <- tcl(canvas_widget, "bbox", "all")
    tcl(canvas_widget, "config", scrollregion = bbox)
  })
  ## adjust "window" width when canvas is resized.
  tkbind(canvas_widget, "<Configure>", function(W) {
    width <- as.numeric(tkwinfo("width", W))
    frame_width <- as.numeric(tkwinfo("width", frame))
    if(frame_width < width)
      tkitemconfigure(canvas_widget, frame_id, width = width)
  })
  return(frame)
}


###################################################
### code chunk number 297: ex-tcltk-scrollable-frame.Rnw:68-73
###################################################
window <- tktoplevel()
tkwm.title(window,"Scrollable frame example")
frame <- ttkframe(window)
tkpack(frame, expand = TRUE, fill = "both")
scroll_frame <- scrollable_frame(frame, 300, 300)


###################################################
### code chunk number 298: ex-tcltk-scrollable-frame.Rnw:82-95
###################################################
font_families <- as.character(tkfont.families())
## skip odd named ones
font_families <- font_families[grepl("^[[:alpha:]]",
                                     font_families)] 
for(i in seq_along(font_families)) {
  font_name <- sprintf("::font::-%s", i)
  try(tkfont.create(font_name, family = font_families[i],
                    size = 14), 
      silent = TRUE)
  l <- ttklabel(scroll_frame, text = font_families[i],
                font = font_name)
  tkpack(l, side = "top", anchor = "w")
}


###################################################
### code chunk number 299: ex-tcltk-sparklines.Rnw:23-30
###################################################
## sparklines
library(tcltk)
library(tseries)
window <- tktoplevel()
tkwm.title(window, "Sparklines example")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 300: ex-tcltk-sparklines.Rnw:35-40
###################################################
mL <- function(label) { # save some typing
  if(is.numeric(label))
    label <- sprintf("%.2f", label)
  ttklabel(frame, text = label, justify = "right") 
}


###################################################
### code chunk number 301: makeHeaderRule
###################################################
tkgrid(mL(""), mL("2000-01-01"), mL("-- until --"), 
       mL("today"), mL("low"), mL("high"))
tkgrid(ttkseparator(frame), row=1, column = 1, columnspan = 5, 
       sticky = "we")


###################################################
### code chunk number 302: ex-tcltk-sparklines.Rnw:56-86
###################################################
add_sparkline <- function(label, symbol = "MSFT") {
  width <- 100; height = 15               # fix width, height
  y <- get.hist.quote(instrument=symbol, start = "2000-01-01",
                      quote = "C", provider = "yahoo", 
                      retclass = "zoo")$Close
  min <- min(y); max <- max(y)
  ##
  start <- y[1]; end <- tail(y,n = 1)
  rng <- range(y)
  ##
  spark_line_canvas <- tkcanvas(frame,
                                width=width, height = height)
  x <- 0:(length(y)-1) * width/length(y)
  if(diff(rng) !=  0) {
    y1 <- (y - rng[1])/diff(rng) * height
    y1 <- height - y1   # adjust to canvas coordinates
  } else {
    y1 <- height/2 + 0 * y
  }
  ## make line with: pathName create line x1 y1... xn yn 
  l <- list(spark_line_canvas, "create","line")
  sapply(seq_along(x), function(i) {
    l[[2*i + 2]] <<- x[i]
    l[[2*i + 3]] <<- y1[i]
  })
  do.call("tcl", l)

  tkgrid(mL(label),mL(start), spark_line_canvas, 
         mL(end), mL(min), mL(max), pady = 2, sticky = "e")
}


###################################################
### code chunk number 303: ex-tcltk-sparklines.Rnw:90-93 (eval = FALSE)
###################################################
## add_sparkline("Microsoft", "MSFT")
## add_sparkline("General Electric", "GE")
## add_sparkline("Starbucks", "SBUX")


###################################################
### code chunk number 304: ex-tcltk-canvas.Rnw:8-14
###################################################
## Canvas example of moving a point. See tkcanvas for much more
library(tcltk)
window <- tktoplevel()
tkwm.title(window, "Move canvas object example")
canvas <- tkcanvas(window, width = 450, height = 300)
tkpack(canvas, expand = TRUE, fill = "both")


###################################################
### code chunk number 305: ex-tcltk-canvas.Rnw:20-25
###################################################
x <- 200; y <- 150; r <- 6
item <- tkcreate(canvas, "oval", x - r, y - r, x + r, y + r,
                 width = 1, outline = "black",
                 fill = "blue")
tkaddtag(canvas, "point", "withtag", item)


###################################################
### code chunk number 306: ex-tcltk-canvas.Rnw:32-36
###################################################
tkitembind(canvas, "point", "<Any-Enter>", function()
           tkitemconfigure(canvas, "current", fill = "red"))
tkitembind(canvas, "point", "<Any-Leave>", function()
           tkitemconfigure(canvas, "current", fill = "blue"))


###################################################
### code chunk number 307: ex-tcltk-canvas.Rnw:42-49
###################################################
last_pos <- numeric(2)            # global to track position
tag_selected <- function(W, x, y) {
  tkaddtag(W,  "selected",  "withtag",  "current")
  tkitemraise(W, "current")
  last_pos <<- as.numeric(c(x, y))
}
tkitembind(canvas, "point", "<Button-1>",  tag_selected)


###################################################
### code chunk number 308: moveSelectedPoint
###################################################
move_selected <- function(W, x, y) {
  pos <- as.numeric(c(x,y))
  tkmove(W, "selected", pos[1] - last_pos[1], 
                        pos[2] - last_pos[2])
  last_pos <<- pos
}
tkbind(canvas, "<B1-Motion>", move_selected)


###################################################
### code chunk number 309: ex-tcltk-canvas.Rnw:70-72
###################################################
tkbind(canvas, "<ButtonRelease-1>", 
       function(W) tkdtag(W, "selected"))


