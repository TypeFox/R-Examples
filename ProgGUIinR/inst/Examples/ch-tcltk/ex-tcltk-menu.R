### R code from vignette source 'ex-tcltk-menu.Rnw'

###################################################
### code chunk number 1: ex-tcltk-menu.Rnw:1-34
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
### code chunk number 2: ex-tcltk-menu.Rnw:40-44
###################################################
library(svMisc)                         # for some helpers
showCmd <- function(cmd) {
  writeLines(captureAll(parseText(cmd)))
}


###################################################
### code chunk number 3: ex-tcltk-menu.Rnw:49-55
###################################################
window <- tktoplevel()
tkwm.title(window, "Simple code editor")
frame <- ttkframe(window, padding = c(3,3,12,12)) 
tkpack(frame, expand = TRUE, fill = "both")
text_buffer <- tktext(frame, undo = TRUE)
addScrollbars(frame, text_buffer)


###################################################
### code chunk number 4: ex-tcltk-menu.Rnw:61-69
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
### code chunk number 5: ex-tcltk-menu.Rnw:75-80
###################################################
tkadd(file_menu, "command", label = "Evaluate buffer",
      command = function() {
        cur_val <- tclvalue(tkget(text_buffer, "1.0", "end"))
        showCmd(cur_val)
      })


###################################################
### code chunk number 6: ex-tcltk-menu.Rnw:84-91
###################################################
tkadd(file_menu, "command", label = "Evaluate selection",
      state = "disabled",
      command =  function() {
        cur_sel <- tclvalue(tkget(text_buffer,
                                  "sel.first", "sel.last"))
        showCmd(cur_sel)
      })


###################################################
### code chunk number 7: addQuit
###################################################
tkadd(file_menu, "separator")
tkadd(file_menu, "command", label = "Quit", 
      command = function() tkdestroy(window))


###################################################
### code chunk number 8: ex-tcltk-menu.Rnw:102-110
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
### code chunk number 9: ex-tcltk-menu.Rnw:116-126
###################################################
tkbind(text_buffer, "<<Selection>>", function(W) {
  hasSelection <- function(W) {
    ranges <- tclvalue(tcl(W, "tag", "ranges", "sel"))
    length(ranges) > 1 || ranges != ""
  }
  ## configure using an index
  sel_state <- ifelse(hasSelection(W), "normal", "disabled")
  print(sel_state)
  tkentryconfigure(file_menu, 2, state = sel_state)
})


###################################################
### code chunk number 10: ex-tcltk-menu.Rnw:130-137
###################################################
tkbind(text_buffer, "<<Modified>>", function(W) {
  ## not really can_undo/can_redo but nothing suitable
  can_undo <- as.logical(tcl(W,"edit", "modified"))
  undo_state <- ifelse(can_undo, "normal", "disabled")
  sapply(c("Undo", "Redo"), function(i)        # match pattern
         tkentryconfigure(edit_menu, i, state = undo_state)) 
})


###################################################
### code chunk number 11: addKeyboardShortcut
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
### code chunk number 12: definePopup
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
### code chunk number 13: addPopup
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


