##' @include guiComponents.R

##' Class for widget to embed a graphic image
setClass("gImage",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' A widget for displaying an image
##'
##' @exports
gimage <- function(
                   filename = "", dirname = "", size = "", handler = NULL,
                   action = NULL, container = NULL, ... ,
                   toolkit=guiToolkit()){
  widget <- .gimage (toolkit,
                     filename=filename, dirname=dirname, size=size,
                     handler=handler, action=action, container=container ,...
                     )
  obj <- new( 'gImage',widget=widget,toolkit=toolkit) 
  return(obj)
}

##' generic for toolkit dispatch
##' @alias gimage
setGeneric( '.gimage' ,
           function(toolkit,
                    filename = "", dirname = "", size = "",
                    handler = NULL,action = NULL, container = NULL, ... )
           standardGeneric( '.gimage' ))


