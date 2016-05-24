##' @include methods.R
NULL





##' A form layout container
##'
##' This convenience container is basically a simpler form of
##' \code{glayout} to be used to layout two columns forms with a label
##' on the left.  The label can be passed in to the \code{add} method
##' of the container as is done with notebook labels
##' @param align alignment of label. Left justify or center
##' balance. Leave as "default" for underlying toolkit default.
##' @param spacing spacing between columns
##' @inheritParams gwidget
##' @export
##' @examples
##' \dontrun{
##' w <- gwindow("gformlayout", visible=FALSE)
##' g <- gvbox(container=w)
##' 
##' flyt <- gformlayout(container=g)
##' gedit("", label="Name:", container=flyt)
##' gedit("", label="Rank:", container=flyt)
##' gedit("", label="Serial No.:", container=flyt)
##' 
##' b <- gbutton("Show me", container=g, handler=function(h,...) {
##' print(svalue(flyt))
##' })
##' 
##' addSpring(g) ## better with Qt, else flyt expands to fill.
##' visible(w) <- TRUE
##' }
gformlayout <- function(
                        align=c("default", "left", "center"),
                        spacing=5,
                        container = NULL, ... ,
                        toolkit=guiToolkit()){

  
  obj <- .gformlayout(toolkit,
                      align=match.arg(align),
                      spacing=spacing,
                      container=container ,...
                      )
  check_return_class(obj, "GFormLayout")
  return(obj)
}


##' .gformlayout generic for toolkit dispatch
##'
##' @export
##' @rdname gformlayout
.gformlayout <- function(toolkit,
                         align="left",
                         spacing=5,
                         container = NULL,
                         ... )
  UseMethod( '.gformlayout' )

##' svalue method
##'
##' The \code{svalue} method for \code{GFormLayout} returns a list of
##' values created by calling \code{svalue} on each child. The
##' returned list is named by the corresponding labels.
##' @inheritParams svalue
##' @export
##' @rdname gformlayout
##' @method svalue GFormLayout
##' @S3method svalue GFormLayout
svalue.GFormLayout <- function(obj, index=NULL, drop=NULL, ...)   NextMethod()



##' svalue assignment method for gformlayout
##'
##' For \code{gformlayout} the \code{svalue} assigment method takes a
##' named list and calls \code{svalue<-} on the children with matching
##' names.
##' @rdname svalue
##' @export
##' @usage \method{svalue}{GFormLayout} (obj, index=NULL, ...) <- value
##' @S3method svalue<- GFormLayout
##' @method svalue<- GFormLayout
"svalue<-.GFormLayout" <- function(obj, index=NULL, ..., value) NextMethod()
