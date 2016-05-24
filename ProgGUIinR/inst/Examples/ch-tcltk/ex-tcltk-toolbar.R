### R code from vignette source 'ex-tcltk-toolbar.Rnw'

###################################################
### code chunk number 1: ex-tcltk-toolbar.Rnw:1-4
###################################################
## Toolbar example in Tk
## tbFrame and Frame are key containers
library(tcltk)


###################################################
### code chunk number 2: ex-tcltk-toolbar.Rnw:18-21
###################################################
window <- tktoplevel(); tkwm.title(window, "Toolbar example")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 3: ex-tcltk-toolbar.Rnw:27-29
###################################################
tool_bar_frame <- ttkframe(frame, padding = 0)
content_frame <- ttkframe(frame, padding = 4)


###################################################
### code chunk number 4: ex-tcltk-toolbar.Rnw:36-44
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
### code chunk number 5: ex-tcltk-toolbar.Rnw:50-52
###################################################
tcl("ttk::style", "configure", "Toolbar.TButton", 
    font = "helvetica 12", padding = 0, borderwidth = 0)


###################################################
### code chunk number 6: ex-tcltk-toolbar.Rnw:58-75
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
### code chunk number 7: addButtons
###################################################
sapply(c("ok", "quit", "cancel"), function(i)
       tkpack(make_icon(tool_bar_frame, i), side = "left"))


###################################################
### code chunk number 8: ex-tcltk-toolbar.Rnw:90-93
###################################################
setState <- function(W, state) tcl(W, "state", state)
tkbind("TButton","<Enter>",function(W) setState(W, "active"))
tkbind("TButton","<Leave>",function(W) setState(W, "!active"))


###################################################
### code chunk number 9: checkStyle
###################################################
function(W) {
  if(as.character(tkcget(W, "-style")) == "Toolbar.TButton")
    cat("... do something for toolbar buttons ...")
}


