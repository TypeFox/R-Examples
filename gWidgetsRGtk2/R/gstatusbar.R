## gtkStatusBar. Use value to push message, value to pop
setClass("gStatusbarRGtk",
         contains="gComponentRGtk",
         representation=representation(label="GtkLabel"),
         prototype=prototype(new("gComponentRGtk"))
         )
## constructor
setMethod(".gstatusbar",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   text="", container=NULL, ...) {

            force(toolkit)
            
            statusbar <- gtkStatusbarNew()
            statusbar$setHasResizeGrip(TRUE)
            sbl <- statusbar[[1]][[1]]
            ## use our own label, not statusbars
            l <- gtkLabel()
            ## l$modifyFont(pangoFontDescriptionFromString("10px"))
            l['xalign'] <- 0.0          # not 0.5
            statusbar[[1]]$remove(statusbar[[1]][[1]])
            statusbar[[1]]$add(l)
            
            ##statusbar$push(statusbar$getContextId("message"), text)

            obj <- as.gWidgetsRGtk2(statusbar)

            svalue(obj) <- text
            
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj,...)
            }
  
            invisible(obj)
          })

as.gWidgetsRGtk2.GtkStatusbar <- function(widget,...) {
  obj <- new("gStatusbarRGtk",block=widget, widget=widget,
             label=widget[[1]][[1]],
             toolkit=guiToolkit("RGtk2"))
  return(obj)
}


### methods

## This pops from stack
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gStatusbarRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            obj@label$getLabel()
            ## obj@widget$Pop(obj@widget$getContextId("message"))
          })

## This pushes to stack
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gStatusbarRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   label <- obj@label
                   if(value == "") value <- " " # need some text to act as strut
                   label$setText(paste(value, collapse="\n"))
                   # obj@widget$Push(obj@widget$getContextId("message"), value)
                   return(obj)
                 })

## push font down to label
setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gStatusbarRGtk"),
                 function(obj, toolkit, ..., value) {
                   .font(obj@label, toolkit, ...) <- value
                   return(obj)
                 })
