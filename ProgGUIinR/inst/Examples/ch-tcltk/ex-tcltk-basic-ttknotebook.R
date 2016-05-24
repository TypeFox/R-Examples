
###################################################
### code chunk number 113: BasicContainers.Rnw:754-758
###################################################
## ttknotebook example
window <- tktoplevel();
tkwm.title(window, "Notebook example")
notebook <- ttknotebook(window)
tkpack(notebook, expand = TRUE, fill = "both")


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
