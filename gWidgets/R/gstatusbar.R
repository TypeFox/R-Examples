##' @include guiComponents.R

##' Class for adding a status bar to main window
setClass("gStatusbar",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' constructor to add a status bar to main window
##'
##' @exports
gstatusbar <- function(
                       text = "", container = NULL, ... ,
                       toolkit=guiToolkit()){
  widget <- .gstatusbar (toolkit,
    text=text, container=container ,...
    )
  obj <- new( 'gStatusbar',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gstatusbar
setGeneric( '.gstatusbar' ,
           function(toolkit,
                    text = "", container = NULL, ... )
           standardGeneric( '.gstatusbar' ))
