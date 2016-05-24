##' @include guiComponents.R

##' a calendar class for date selection
setClass("gCalendar",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' A constructor for a date selection widget
##'
##' @param text initial text
##' @param format Date format
##' @param handler handler called when changed
##' @param action passed to handler
##' @param container parent container
##' @param ... passed to \code{add} method of parent
##' @param toolkit toolkit
##' @return Returns an object of class \code{gCalendar} for which the following methods are overridden:
##' \enumerate{
##' \item \code{svalue} get the date
##' 
##' \item \code{svalue<-} set the date
##' }
##' The change handler is inherited from \code{\link{gedit}}
##' @export
gcalendar <- function(
                      text = "", format = "%Y-%m-%d", 
                      handler = NULL, action=NULL, container = NULL,...,
                      toolkit=guiToolkit()){
  widget <- .gcalendar (toolkit,
                        text=text, format=format, handler=handler,action=action,
                        container=container , ...
                        )
  obj <- new( 'gCalendar',widget=widget,toolkit=toolkit) 
 return(obj)
}


##' generic for toolkit dispatch
##' @alias gcalendar
setGeneric( '.gcalendar' ,
           function(toolkit,
                    text = "", format = "%Y-%m-%d", 
                    handler=NULL, action=NULL, container = NULL,
                    ... )
           standardGeneric( '.gcalendar' ))



##' svalue method for gcalendar
##'
##' Main property of calendar is the date
##' @param obj object
##' @param index ignored
##' @param drop ignored
##' @return character. The date, after formatting.
##' @exports
setMethod("svalue", signature(obj="gCalendar"),
          function(obj, index=NULL, drop=NULL, ... ) {
            .svalue(obj@widget, obj@toolkit, ...,index=index, drop=drop)            
          })




##' set date
##'
##' @param obj
##' @param index ignored
##' @param ... ignored
##' @param value character. Should be in format for calendar widget
##' @return void
##' @exports
setReplaceMethod("svalue", signature(obj="gCalendar"),
          function(obj, index=NULL, ...,value) {
            .svalue(obj@widget, obj@toolkit, index=index, ...) <- value
            return(obj)
          })

