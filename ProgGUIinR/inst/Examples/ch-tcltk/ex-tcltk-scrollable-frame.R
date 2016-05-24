### R code from vignette source 'ex-tcltk-scrollable-frame.Rnw'

###################################################
### code chunk number 1: ex-tcltk-scrollable-frame.Rnw:1-20
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
### code chunk number 2: ex-tcltk-scrollable-frame.Rnw:39-64
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
### code chunk number 3: ex-tcltk-scrollable-frame.Rnw:68-73
###################################################
window <- tktoplevel()
tkwm.title(window,"Scrollable frame example")
frame <- ttkframe(window)
tkpack(frame, expand = TRUE, fill = "both")
scroll_frame <- scrollable_frame(frame, 300, 300)


###################################################
### code chunk number 4: ex-tcltk-scrollable-frame.Rnw:82-95
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


