##' @include guiContainer.R

##' Class for tabbed notebook container
setClass("gNotebook",
         contains="guiContainer",
         prototype=prototype(new("guiContainer"))
         )

##' Constructor for a tabbed notebook container
##'
##' @export
gnotebook <- function(
                      tab.pos = 3, closebuttons = FALSE, dontCloseThese = NULL,
                      container = NULL, ... ,
                      toolkit=guiToolkit()){
  widget <- .gnotebook (toolkit,
                        tab.pos=tab.pos, closebuttons=closebuttons, dontCloseThese=dontCloseThese,
                        container=container ,...
                        )
  obj <- new( 'gNotebook', widget=widget,toolkit=toolkit)
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gnotebook
setGeneric( '.gnotebook' ,
           function(toolkit,
                    tab.pos = 3, closebuttons = FALSE, dontCloseThese = NULL,
                    container = NULL, ... )
           standardGeneric( '.gnotebook' ))
