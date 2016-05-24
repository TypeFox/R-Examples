## Could make spinbutton slider, subclass as methods are identical
setClass("gSpinbuttonRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

setMethod(".gspinbutton",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   from=0,to=10,by=1,value=from,digits=0,
                   handler=NULL, action=NULL,
                   container=NULL, ...) {

            force(toolkit)

            ## fix digits if user forgot
            if(digits == 0 &&  as.logical((by %% 1))) # FALSE if integer o/w T
              digits = abs(floor(log(by,10)))
             
            
            adjustment = gtkAdjustmentNew(value=value, lower=from,
              upper=to,step.incr=by)
            spin = gtkSpinButtonNew(adjustment,climb.rate=0.6, digits=digits)

            obj <- as.gWidgetsRGtk2(spin) 

            svalue(obj) <- value                  # wasn't working as desired
  

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

as.gWidgetsRGtk2.GtkSpinButton <- function(widget,...) {
  obj <- new("gSpinbuttonRGtk", block=widget, widget=widget,
             toolkit=guiToolkit("RGtk2"))
  return(obj)
}

### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gSpinbuttonRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            obj@widget$GetValue()
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gSpinbuttonRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   obj@widget$SetValue(value)
                   return(obj)
                 })

## Method to replace values of sping button
setReplaceMethod("[",
                 signature(x="gSpinbuttonRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gSpinbuttonRGtk"),
          function(x, toolkit, i, j, ..., value) {
            obj <- x
            widget <- getWidget(obj)

            ## check that value is a regular sequence
            if(length(value) <=1) {
              warning("Can only assign a vector with equal steps, as produced by seq")
              return(obj)
            }
            if(length(value) > 2 &&
               !all.equal(diff(diff(value)), rep(0, length(value) - 2))) {
              warning("Can only assign a vector with equal steps, as produced by seq")
              return(obj)
            }
            ## get current value, increment
            curValue <- svalue(obj)
            inc <- head(diff(value), n=1)

            widget$setRange(min(value), max(value))
            widget$setIncrements(inc, inc) # button 1, button 2
            widget$setValue(curValue)


            ## all done
            return(obj)
          })




### handlers
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gSpinbuttonRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandler(obj,"value-changed",handler, action, ...)
          })

