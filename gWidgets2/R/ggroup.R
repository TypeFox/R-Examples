##' @include methods.R
NULL

##' Basic box container
##'
##' @param horizontal logical. If TRUE, left to right layout, otherwise top to bottom
##' @param spacing spacing aroud widget 
##' @param use.scrollwindow logical. Either \code{TRUE},
##' \code{"TRUE"}, \code{FALSE}, \code{"FALSE"}, \code{"y"}, or
##' \code{"x"}. For all toolkits a non-FALSE value will place the
##' child components into a scrollable container. For some toolkits
##' this will only be in the direction of packing. If the toolkit
##' allows it (RGtk2), then values of \code{"x"} or \code{"y"} can be
##' used to override the default scrolling directions. A box container with
##' scrollwindows should have it size set either directly or through
##' packing with \code{expand=TRUE} as its size request will not
##' reflect the size of its child components.
##' 
##' @inheritParams gwidget
##' @return a GGroup instance.
##' @export
##' @rdname ggroup
##' @seealso \code{\link{gframe}} and \code{\link{gexpandgroup}}
##' @example inst/examples/ex-boxcontainers.R
ggroup <- function(horizontal=TRUE, spacing=5, use.scrollwindow=FALSE, container=NULL, ..., toolkit=guiToolkit()) {

  if(is.character(toolkit))
    toolkit <- guiToolkit(toolkit)
  
  obj <- .ggroup(toolkit, horizontal, spacing=spacing, use.scrollwindow=use.scrollwindow, container, ...)

  check_return_class(obj, "GGroup")
  obj   
  
}

##' S3 generic whose methods are implemented in the toolkit packages
##'
##' @rdname ggroup
##' @export
.ggroup <- function(toolkit, horizontal=TRUE, spacing=5, use.scrollwindow=FALSE, container=NULL, ...) UseMethod(".ggroup")

##' \code{svalue<-} method for a ggroup
##'
##' The \code{svalue} method refers to the main property of the box
##' container, its spacing. There are generally two types of spacing:
##' padding around border of the box and spacing between each child
##' that is packed in. The spacing here is the between-child-component spacing.
##' The reference class method \code{set_borderwidth} can be used for the other.
##'
##' Child components are typically added to a box container through
##' the child components constructor. The argument \code{expand},
##' \code{fill}, and \code{anchor} determine how the child is
##' positioned within the container.
##' @usage \method{svalue}{GGroup} (obj, index=TRUE, ...) <- value
##' @param obj \code{GGroup} object
##' @param index ignored
##' @param value value (in pixels) for between child spacing
##' @export
##' @rdname ggroup
##' @method svalue<- GGroup
##' @S3method svalue<- GGroup
"svalue<-.GGroup" <- function(obj, index=TRUE,  ..., value) {
  NextMethod()
}

##' Convenience constructor for vertical ggroup
##'
##' Avoids need to type \code{horizontal=FALSE}
##' @inheritParams ggroup
##' @return a GGroup instance with vertical packing.
##' @export
##' @rdname ggroup
gvbox <- function(spacing=5, use.scrollwindow=FALSE, container=NULL, ..., toolkit=guiToolkit())
  ggroup(horizontal=FALSE, spacing=spacing, use.scrollwindow=use.scrollwindow, container=container, ..., toolkit=toolkit)
