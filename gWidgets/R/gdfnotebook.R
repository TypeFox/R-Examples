##' @include guiComponents.R

##' class to hold a notebook of data frame editors
setClass("gDfNotebook",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' A notebook container for many \code{gdf} instances
##'
##' @exports
##' @param items data frame for initial page
##' @param container parent container
##' @param ... passed to \code{add} method of parent container
##' @param toolkit toolkit
gdfnotebook <- function(
                        items = NULL, container = NULL, ... ,
                        toolkit=guiToolkit()){
  widget <- .gdfnotebook (toolkit,
                          items=items, container=container ,...
                          )
  obj <- new( 'gDfNotebook',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gdfnotebook
setGeneric( '.gdfnotebook' ,
           function(toolkit,
                    items = NULL, container = NULL, ... )
           standardGeneric( '.gdfnotebook' ))
