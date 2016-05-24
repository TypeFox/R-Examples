##' @include gframe.R

##' Class for a box container with disclosure trigger
setClass("gExpandGroup",
         contains="gFrame",
         prototype=prototype(new("gFrame"))
         )

##' Constructor of box container widget with disclosure trigger and label
##'
##' @export
##' @param text Label text
##' @param markup logical. Does text have markup. (not too many)
##' @param horizontal horizontal (\code{TRUE}) or vertical packing.
##' @param handler handler called when toggled
##' @param action passed to handler
##' @param container parent container
##' @param ... passed to parent's \code{add} method
##' @param toolkit toolkit
##' @return An object of class \code{gExpandgroup}. This (basically)
##' inherits from \code{gFrame} its methods and overrides:
##'
##' \enumerate{
##'
##' \item \code{visible<-} Logical. To specify if widget is open (\code{TRUE}) or closed.
##'
##' }
##'
##' 
gexpandgroup <- function(
                         text = "", markup = FALSE, horizontal=TRUE,
                         handler = NULL, action = NULL,
                         container = NULL, ... ,
                         toolkit=guiToolkit()){
  widget <- .gexpandgroup (toolkit,
                           text=text, markup=markup, horizontal=horizontal,
                           handler=handler, action=action, container=container ,...
                           )
  obj <- new( 'gExpandGroup',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gexpandgroup
setGeneric( '.gexpandgroup' ,
           function(toolkit,
                    text = "", markup = FALSE,horizontal=TRUE,
                    handler = NULL, action = NULL,
                    container = NULL, ... )
standardGeneric( '.gexpandgroup' ))
