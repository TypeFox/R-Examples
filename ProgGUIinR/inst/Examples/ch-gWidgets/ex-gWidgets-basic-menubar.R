###################################################
### code chunk number 133: Controls.Rnw:1424-1441
###################################################
stub <- function(h,...) gmessage("called handler", 
                                 parent = window)
action_list = list(
  new = gaction(label = "new", icon = "new", 
    handler = stub, parent = window),
  open = gaction(label = "open", icon = "open", 
    handler = stub, parent = window),
  save = gaction(label = "save", icon = "save", 
    handler = stub, parent = window),
  save.as = gaction(label = "save as...", icon = "save as...", 
    handler = stub, parent = window),
  quit = gaction(label = "quit", icon = "quit", 
    handler = function(...) dispose(window), parent = window),
  cut = gaction(label = "cut", icon = "cut", 
    handler = stub, parent = window)
  )



###################################################
### code chunk number 134: Controls.Rnw:1445-1451
###################################################
window <- gwindow("gtoolbar example")
tool_bar_list<- c(action_list[c("new","save")], 
                 sep = gseparator(), 
                 action_list["quit"])
tool_bar <- gtoolbar(tool_bar_list, cont = window)
gtext("Lorem ipsum ...", cont = window)


###################################################
### code chunk number 135: Controls.Rnw:1478-1491
###################################################
menu_bar_list <- list(file = list(
             new = action_list$new,
             open = action_list$open,
             save = action_list$save,
             "save as..." = action_list$save.as,
             sep = gseparator(),
             quit = action_list$quit
             ),
           edit = list(
             cut = action_list$cut
             )
           )




###################################################
### code chunk number 136: Controls.Rnw:1494-1498
###################################################
window <- gwindow("Menu bar example")
menu_bar <- gmenu(menu_bar_list, cont = window)
tool_bar <- gtoolbar(tool_bar_list, cont = window)
txt_widget <- gtext("", cont = window, expand = TRUE)


###################################################
### code chunk number 137: Controls.Rnw:1527-1535
###################################################
no_changes <- c("save","save.as","cut")
keyhandler <- function(...) {
  for(i in no_changes)
    enabled(action_list[[i]]) <- 
      (nchar(svalue(txt_widget)) > 0)
}
addHandlerKeystroke(txt_widget, handler = keyhandler)
keyhandler()
