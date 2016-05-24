
###################################################
### code chunk number 43: tabbed_notebook
###################################################
window <- gwindow("gnotebook example")
notebook <- gnotebook(cont = window)


###################################################
### code chunk number 44: addPage
###################################################
add_a_page <- function(file_name) {
  f <- system.file(file_name, package = "ProgGUIinR")
  gtext(paste(readLines(f), collapse="\n"), cont = notebook, label = file_name)
}
add_a_page("DESCRIPTION")


###################################################
### code chunk number 45: helpPage
###################################################
lyt <- glayout(cont = notebook, horizontal = FALSE, 
               label = "Help")
lyt[1,1] <- gimage("help", dir = "stock", cont = lyt)
lyt[1,2] <- glabel(paste("To add a page:",
           "Click on a file in the left pane, and its contents",
           "are displayed in a notebook page.", sep = "\n"), 
           cont = lyt)


###################################################
### code chunk number 46: set_notebook_page
###################################################
svalue(notebook) <- length(notebook)

