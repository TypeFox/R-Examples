##' @include guiComponents.R

##' Class for display of line in a layout
setClass("gSeparator",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' constructor providing a widget for displaying a line in a GUI
##'
##' @exports
gseparator <- function(
                       horizontal = TRUE, container = NULL, ... ,
                       toolkit=guiToolkit()){
  widget <- .gseparator (toolkit,
                         horizontal=horizontal, container=container ,...
                         )
  obj <- new( 'gSeparator',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gseparator
setGeneric( '.gseparator' ,
           function(toolkit,
                    horizontal = TRUE, container = NULL, ... )
           standardGeneric( '.gseparator' ))
