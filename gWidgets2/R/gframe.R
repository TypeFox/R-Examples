##' @include ggroup.R
NULL

##' Constructor for framed box container with label
##'
##' The framed box container inherits from \code{ggroup}. The main
##' addition is a label, which is accessed via the \code{name} method.
##' @param text frame label
##' @param markup does label use markup (toolkit specific)
##' @param pos position of label: 0=left, 1=right, some toolkit allow values in between
##' @param horizontal logical. If TRUE, left to right layout, otherwise top to bottom
##' @param spacing spacing aroud widget 
##' @param container parent container
##' @param ... passed through
##' @param toolkit toolkit
##' @seealso \code{\link{ggroup}} and \code{\link{gexpandgroup}}
##' @note to include a scrollwindow, place a \code{ggroup} within this window.
##' @export
##' @rdname gframe
##' @examples
##' \dontrun{
##' w <- gwindow("gformlayout", visible=FALSE)
##' f <- gframe("frame", horizontal=FALSE, container=w)
##' glabel("Lorem ipsum dolor sit amet, \nconsectetur adipiscing elit.", container=f)
##' gbutton("change name", container=f, handler=function(h,...) {
##'   names(f) <- "new name"
##' })
##' visible(w) <- TRUE
##' }
gframe <- function(
                   text = "", markup=FALSE, pos = 0, horizontal=TRUE, spacing=5, container = NULL,
                   ... ,
                   toolkit=guiToolkit()){
  obj <- .gframe (toolkit,
           text,  markup, pos, horizontal, spacing, container,
           ...
           )

  check_return_class(obj, "GFrame")
  obj   
  
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gframe
.gframe <- function(toolkit,
                    text = "", markup = FALSE, pos = 0, horizontal=TRUE, spacing=5,
                    container = NULL,      ... )
           UseMethod( '.gframe' )

##' set names for frame
##'
##' @export
##' @usage \method{names}{GFrame} (x) <- value
##' @rdname gWidgets2-S3methods
##' @method names<- GFrame
##' @S3method names<- GFrame
"names<-.GFrame" <- function(x, value) NextMethod()
