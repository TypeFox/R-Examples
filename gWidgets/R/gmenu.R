##' @include guiComponents.R

##' single line text edit class
setClass("gMenu",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' menu constructor, main and popup
##'
##' @exports
gmenu <- function(
  menulist, popup = FALSE, action = NULL, container = NULL,      ... ,
  toolkit=guiToolkit()){
  widget <- .gmenu (toolkit,
                    menulist=menulist, popup=popup, action=action, container=container ,...
                    )
  obj <- new( 'guiComponent',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gmenu
setGeneric( '.gmenu' ,
           function(toolkit,
                    menulist, popup = FALSE, action = NULL, container = NULL,
                    ... )
           standardGeneric( '.gmenu' ))
