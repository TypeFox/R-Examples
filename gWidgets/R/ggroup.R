##' @include guiContainer.R

##' Base class for box containers. Subclasses are gFrame, gExpandGroup
setClass("gGroup",
         contains="guiContainer",
         prototype=prototype(new("guiContainer"))
         )

##' Constructor for horizontal or vertical box container
##'
##' @param horizontal layout children left to right (\code{TRUE}) or top to bottom (\code{FALSE})
##' @param spacing between widget spacing. (No API to set margin size)
##' @param use.scrollwindow if \code{TRUE} uses scrollbars when
##' requested size for children larger than display size of the widget
##' @param container parent container
##' @param ... ignored
##' @param toolkit underlying toolkit, unusual to specify
##' @export
##' @seealso \code{\link{gframe}}, \code{\link{gexpandgroup}}
ggroup <- function(
                   horizontal = TRUE, spacing = 5, use.scrollwindow = FALSE, container = NULL, ... ,
                   toolkit=guiToolkit()){
  widget <- .ggroup (toolkit,
                     horizontal=horizontal, spacing=spacing,
                     use.scrollwindow = use.scrollwindow, 
                     container=container,...
                     )
  obj <- new( 'gGroup',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' 
##' @alias ggroup
setGeneric( '.ggroup' ,
           function(toolkit,
                    horizontal = TRUE, spacing = 5,  use.scrollwindow = FALSE, container = NULL,  ... )
           standardGeneric( '.ggroup' ))



################## Methods ###############################


##' Return spacing between widgets in a box container
##'
##' Main property of a box container is the spacing between widgets
##' @param obj object
##' @param index ignored
##' @param drop ignored
##' @param ... ignored
##' @return integer. between widget spacing in pixes
setMethod("svalue",signature(obj="gGroup"),
          function(obj, index=NULL, drop=NULL, ... ) {
            .svalue(obj@widget, obj@toolkit, ...,index=index, drop=drop)            
          })


##' set between widget spacing for box containers
##'
##' @param obj
##' @param index ignormed
##' @param ... ignored
##' @param value integer, spacing values
setReplaceMethod("svalue",signature(obj="gGroup"),
          function(obj, index=NULL, ...,value) {
            .svalue(obj@widget, obj@toolkit, index=index, ...) <- value
            return(obj)
          })


################## addSpace ################################
##' addSpace generic for box containers
##'
##' @export
setGeneric("addSpace",function(obj,value, ...) standardGeneric("addSpace"))

##' Add space between previous child and next one
##'
##' @param obj object
##' @param value integer, pixel size
##' @param ... ignored
##' @return void
setMethod("addSpace",signature(obj="gGroup"),
          function(obj, value, ...) {
            toolkit = obj@toolkit
            .addSpace(obj@widget,toolkit,value,...)
          })

##' dispatch with toolkit
##' @alias addSpace
setGeneric(".addSpace",function(obj,toolkit,value,...) standardGeneric(".addSpace"))


##################################################
##' addSpring generic
setGeneric("addSpring",function(obj,...) standardGeneric("addSpring"))

##' addSpring between children for box containers
##'
##' A spring pushes the two children to the sides (or top and bottom) of the box
##' @param obj object
##' @param ... ignored
##' @return void
setMethod("addSpring",signature(obj="gGroup"),
          function(obj, ...) {
            toolkit = obj@toolkit
            .addSpring(obj@widget, toolkit,...)
          })

##' dispatch with toolkit
##' @alias addSpring
setGeneric(".addSpring",function(obj, toolkit,...) standardGeneric(".addSpring"))


