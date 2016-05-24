##' @include guiContainer.R

##' Class for a two-paned container
setClass("gPanedGroup",
         contains="guiContainer",
         prototype=prototype(new("guiContainer"))
         )

##' constructor for a two-paned container
##'
##' @export
gpanedgroup <- function(
  widget1=NULL, widget2=NULL, horizontal = TRUE, container = NULL , ...,
  toolkit=guiToolkit()){
  widget <- .gpanedgroup (toolkit,
                          widget1, widget2, horizontal=horizontal,
                          container=container, ...
                          )
  obj <- new( 'gPanedGroup',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gpanedgroup
setGeneric( '.gpanedgroup' ,
           function(toolkit,
                    widget1, widget2, horizontal = TRUE, container = NULL, ... )
           standardGeneric( '.gpanedgroup' ))
