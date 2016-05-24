
###################################################
### code chunk number 74: Containers.Rnw:655-673
###################################################
gtkNotebookInsertPageWithCloseButton <- 
  function(object, child, label.text="", position=-1) {
    icon <- gtkImage(pixbuf = 
      object$renderIcon("gtk-close", "button", size = "menu"))
    closeButton <- gtkButton()
    closeButton$setImage(icon)
    closeButton$setRelief("none")
    ##
    label <- gtkHBox()
    label$packStart(gtkLabel(label.text))
    label$packEnd(closeButton)
    ##
    gSignalConnect(closeButton, "clicked", function(button) {
      index <- object$pageNum(child)
      object$removePage(index)
    })
    object$insertPage(child, label, position)
  }


###################################################
### code chunk number 75: Containers.Rnw:678-684
###################################################
window <- gtkWindow()
notebook <- gtkNotebook(); window$add(notebook)
notebook$insertPageWithCloseButton(gtkButton("hello"), 
                                   label.text = "page 1")
notebook$insertPageWithCloseButton(gtkButton("world"), 
                                   label.text = "page 2")

