##' @include guiComponents.R

##' Class for graphic device widget
setClass("gGraphics",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' Constructor for an embeddable graphics device
##'
##' @exports
ggraphics <- function(
                      width = dpi * 6, height = dpi * 6, dpi = 75, ps = 12,      container = NULL, ... ,
                      toolkit=guiToolkit()){
  widget <- .ggraphics (toolkit,
                        width=width, height=height, dpi=dpi, ps=ps, container=container ,...
                        )
  obj <- new( 'gGraphics',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias ggraphics
setGeneric( '.ggraphics' ,
           function(toolkit,
                    width = dpi * 6, height = dpi * 6, dpi = 75, ps = 12,
                    container = NULL, ... )
           standardGeneric( '.ggraphics' ))



##' Class for notebook to hold multiple ggraphics widgets
setClass("gGraphicsNotebook",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' constructor for notebook to hold multiple graphics devices
##'
##' @export
ggraphicsnotebook <- function(
                              width = dpi * 6, height = dpi * 6, dpi = 75, container = NULL,      ... ,
                              toolkit=guiToolkit()){
  widget <- .ggraphicsnotebook (toolkit,
                                width=width, height=height, dpi=dpi, container=container ,...
                                )
  obj <- new( 'gGraphicsNotebook',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias ggraphicsnotebook
setGeneric( '.ggraphicsnotebook' ,
           function(toolkit,
                    width = dpi * 6, height = dpi * 6, dpi = 75,
                    container = NULL,      ... )
           standardGeneric( '.ggraphicsnotebook' ))
