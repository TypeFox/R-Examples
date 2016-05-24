## So much is identical here to gspinbutton, we should make a class to derive these from -- another day.


setClass("gSliderRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

setMethod(".gslider",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   from=0, to=100, by = 1,
                   value=from,
                   horizontal=TRUE,
                   handler=NULL, action=NULL,
                   container=NULL, ...) {

            force(toolkit)

            if(length(from) == 1)
              x <- seq(from, to, by)
            else
              x <- from

            x <- sort(unique(x))
  
            if (horizontal)
              widget <- gtkHScaleNewWithRange(1L, length(x), 1L)
            else
              widget <- gtkVScaleNewWithRange(1L, length(x), 1L)


            
            obj <- as.gWidgetsRGtk2(widget)
            
            ## obj <- new("gSliderRGtk",block=align, widget=widget,
            ##            toolkit=guiToolkit("RGtk2"))
            
            tag(obj, "..byIndexValues") <- x
            tag(obj, "default_fill") <- ifelse(horizontal, "x", "y")
            
            svalue(obj) <- value[1]
            
            gSignalConnect(widget, "format-value", function(widget, value, ...) {
              format(tag(obj, "..byIndexValues")[as.integer(value)], digits=3)
            })
            
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj,...)
            }
            
            if (!is.null(handler))  {
              id = addhandlerchanged(obj, handler, action)
            }
            
            invisible(obj)
          })

##' coercoe gtkwidget into scale widget so that methods can work
as.gWidgetsRGtk2.GtkHScale <-  function(widget, ...) {
  asgWidgetsRGtk2.SCALE(widget, yscale=0, ...)
}
as.gWidgetsRGtk2.GtkVScale <- function(widget, ...) {
  asgWidgetsRGtk2.SCALE(widget, xscale=0, ...)
}

asgWidgetsRGtk2.SCALE <- function(widget,xscale=1, yscale=1, ...) {
    if(is.null(widget$parent)) {
      align <- gtkAlignmentNew(xscale=xscale, yscale=yscale)
      align$add(widget)
      obj <- new("gSliderRGtk",block=align, widget=widget,
                 toolkit=guiToolkit("RGtk2"))
    } else {
      obj <- new("gSliderRGtk",block=widget, widget=widget,
                 toolkit=guiToolkit("RGtk2"))
    }
    return(obj)
  }

### methods

setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gSliderRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            ind <- obj@widget$getValue()
            if(!is.null(index) && index)
              return(ind)
            else
              return(tag(obj, "..byIndexValues")[ind])
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gSliderRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   if(is.null(index) || !index) {
                     ## value is a value, must match
                     value <- as.character(match(value, tag(obj, "..byIndexValues")))
                   }
                   getWidget(obj)$setValue(value)

                   ## update label?
                   return(obj)
                 })

##' return values
##' @param i, j, drop ignored
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gSliderRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            tag(x, "..byIndexValues")
          })

##' non-essential method to dispatch done to leftBracket
setReplaceMethod("[",
                 signature(x="gSliderRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })




setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gSliderRGtk"),
          function(x, toolkit, i, j, ..., value) {
            obj <- x
            widget <- getWidget(obj)
            curValue <- svalue(obj)

            value <- sort(unique(value))
            tag(obj, "..byIndexValues") <- value
            
            widget$setRange(1, length(value))
            widget$setIncrements(1L, 1L) # button 1, button 2

            svalue(obj) <- curValue


            ## all done
            return(obj)
          })



### handlers
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gSliderRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandler(obj, "value-changed", handler, action,...)
          })

