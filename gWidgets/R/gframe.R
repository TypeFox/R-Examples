##' @include ggroup.R

##' Class for framed box container with label
setClass("gFrame",
         contains="gGroup",
         prototype=prototype(new("gGroup"))
         )

##' Constructor for framed box container with label
##'
##' @export
gframe <- function(
                   text = "", markup = FALSE, pos = 0, horizontal=TRUE, container = NULL,
                   ... ,
                   toolkit=guiToolkit()){
  widget <- .gframe (toolkit,
                     text=text, markup=markup, pos=pos, horizontal=horizontal, container=container ,
                     ...
                     )
  obj <- new( 'gFrame',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gframe
setGeneric( '.gframe' ,
           function(toolkit,
                    text = "", markup = FALSE, pos = 0, horizontal=TRUE,
                    container = NULL,      ... )
           standardGeneric( '.gframe' ))
