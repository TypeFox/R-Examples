setClass("gCheckboxRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

## constructor
setMethod(".gcheckbox",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   text, checked=FALSE,
                   use.togglebutton=FALSE,
                   handler=NULL, action=NULL,
                   container=NULL,...) {

            force(toolkit)

            if(missing(text))
              text <- ""

            if(use.togglebutton)
              return(gtogglebutton(text, checked, handler, action, container, ...))
            
            
            check <- gtkCheckButtonNewWithLabel(text)
            check$SetActive(checked)

            obj <- as.gWidgetsRGtk2(check)

            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE, toolkit=toolkit)
              add(container, obj,...)
            }
            
            
            if (!is.null(handler)) {
              id = addhandler(obj, "toggled",handler, action=action)
            }
  
            invisible(obj)
          })

as.gWidgetsRGtk2.GtkCheckButton <- function(widget,...) {
  parent <- widget$parent
  if(is.null(parent)) {
    parent <- gtkAlignmentNew(0,0,0,0)
    parent$add(widget)
  }
  obj = new("gCheckboxRGtk",block=parent, widget=widget,
    toolkit=guiToolkit("RGtk2"))

  return(obj)
}
### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            obj@widget$getActive()
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   obj@widget$setActive(value)
                   return(obj)
                 })

## [
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            x@widget[[1]]$GetText()
          })
            
setMethod("[",
          signature(x="gCheckboxRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gCheckboxRGtk"),
          function(x, toolkit, i, j, ..., value) {
            x@widget[[1]]$SetText(value)
            return(x)
          })

setReplaceMethod("[",
                 signature(x="gCheckboxRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

## handlers
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandler(obj, "toggled", handler, action=action,...)
          })

setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCheckboxRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerclicked(obj, toolkit, handler, action=action,...)
          })


##################################################
##' class to provide a toggle button alternative to a checkbox. The toggle button
##' is very similar
setClass("gToggleButtonRGtk",
         contains="gCheckboxRGtk"
         )

##' Provides a toggle button alternative to a check box.
##' 
##' constructor, not a method as called internally
gtogglebutton <- function(text, checked=FALSE, handler=NULL, action=NULL,
                          container=NULL, ...) {

  widget <- gtkToggleButton()
  
  obj <- new("gToggleButtonRGtk", 
             block=widget, widget=widget,
             toolkit=guiToolkit("RGtk2"))

  if(!missing(text))
    obj[] <- text

  svalue(obj) <- checked
  
  if(!is.null(handler))
    addHandlerChanged(obj, handler=handler, action=action)

  if(!is.null(container)) {
    if(is.logical(container) && container)
      container <- gwindow()

    add(container, obj, ...)

  }

  return(obj)
}

## method to set text
setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gToggleButtonRGtk"),
          function(x, toolkit, i, j, ..., value) {
            tb <- getWidget(x)
            
            icons <- getStockIcons()
            if(value %in% names(icons)) {
              tb['use-stock'] <- TRUE
              tb['label'] <- icons[[value]]
            } else {
              tb['use-stock'] <- FALSE
              tb['label'] <- value
            }

            return(x)
          })
