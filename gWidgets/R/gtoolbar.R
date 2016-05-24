##' @include guiComponents.R

##' A toolbar class
setClass("gToolbar",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' A toolbar constructor
##'
##' @exports
gtoolbar <- function(
                     toolbarlist, style = c("both", "icons", "text", "both-horiz"),
                     action = NULL, container = NULL, ... ,
                     toolkit=guiToolkit()){
  widget <- .gtoolbar (toolkit,
                       toolbarlist=toolbarlist, style=style, action=action, container=container ,...
                       )
  obj <- new( 'gToolbar',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gtoolbar
setGeneric( '.gtoolbar' ,
           function(toolkit,
                    toolbarlist, style = c("both", "icons", "text", "both-horiz"),
                    action = NULL, container = NULL, ... )
           standardGeneric( '.gtoolbar' ))

