##################################################
## add a separator to a container. Needs the container

setClass("gSeparatorRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

## should this return object?
setMethod(".gseparator",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   horizontal = TRUE, container = NULL, ...) {

            force(toolkit)
            
            if(horizontal) {
              separator = gtkHSeparatorNew()
            } else {
              separator = gtkVSeparatorNew()
            }

            obj <- as.gWidgetsRGtk2(separator)

            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj,...)
            }

            invisible(obj)
            
          })

as.gWidgetsRGtk2.GtkHSeparator <- as.gWidgetsRGtk2.GtkVSeparator <- function(widget,...){
  obj <- new("gSeparatorRGtk", block=widget, widget=widget,
             toolkit=guiToolkit("RGtk2"))
  return(obj)
}

.isgSeparator <- function(obj) {
  (is(obj,"guiComponent") && is(obj@widget,"gSeparatorRGtk") ) ||
    is(obj,"gSeparatorRGtk")
}
