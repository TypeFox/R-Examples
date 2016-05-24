### R code from vignette source 'ex-tcltk-tree.Rnw'

###################################################
### code chunk number 1: ex-tcltk-tree.Rnw:3-24
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
### code chunk number 2: ex-tcltk-tree.Rnw:45-50
###################################################
library(XML)
file_name <- "http://www.omegahat.org/RSXML/shortIntro.html"
doc <- htmlTreeParse(file_name, useInternalNodes = TRUE, 
                     error = function(...) {})
root <- xmlRoot(doc)


###################################################
### code chunk number 3: notShown
###################################################
window <- tktoplevel()
tkwm.title(window, "Treeview example with XML")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 4: ex-tcltk-tree.Rnw:61-64
###################################################
treeview <- ttktreeview(frame, displaycolumns = "#all", 
                        columns = 1)
addScrollbars(frame, treeview)                    


###################################################
### code chunk number 5: columnConfiguration
###################################################
tcl(treeview, "heading", "#0", text = "Name")
tcl(treeview, "column",  "#0", minwidth = 20)
tcl(treeview, "heading",  1,   text = "value")
tcl(treeview, "column",   1,   minwidth = 20)


###################################################
### code chunk number 6: ex-tcltk-tree.Rnw:75-95
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
### code chunk number 7: ex-tcltk-tree.Rnw:104-118
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
### code chunk number 8: ex-tcltk-tree.Rnw:131-133
###################################################
.selected_id <- ""                               # globals
.dragging <- FALSE


###################################################
### code chunk number 9: ex-tcltk-tree.Rnw:137-140
###################################################
tkbind(treeview, "<Button-1>", function(W,x,y) {
  .selected_id <<- as.character(tcl(W, "identify","row",x, y))
})  


###################################################
### code chunk number 10: ex-tcltk-tree.Rnw:145-149
###################################################
tkbind(treeview, "<B1-Motion>", function(W,x, y, X, Y) {
  tkconfigure(W, cursor = "diamond_cross")
  .dragging <<-TRUE
})


###################################################
### code chunk number 11: ex-tcltk-tree.Rnw:159-169
###################################################
tkbind(treeview, "<ButtonRelease-1>", function(W,x, y, X, Y) {
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
### code chunk number 12: walkTreeReturnAList
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



