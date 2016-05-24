##' @include guiComponents.R

##' Class for display of HTML text
setClass("gHtml",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' constructor of widget to display HTML text
##'
##' @export
ghtml <- function(
                  x, handler = NULL, 
                  action = NULL, container = NULL, 
                  ..., toolkit=guiToolkit()) {
  widget <- .ghtml(toolkit,
                   x,
                   handler = handler, 
                   action = action, container = container, 
                   ...)
  obj <- new("gHtml",widget=widget,toolkit=toolkit)
  return(obj)
}

##' generic for toolkit dispatch
##' @alias ghtml
setGeneric(".ghtml",function(toolkit,
                             x, handler = NULL, 
                             action = NULL, container = NULL, 
                             ...) standardGeneric(".ghtml"))



