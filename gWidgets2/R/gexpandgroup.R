##' @include gframe.R
NULL

##' Constructor of box container widget with disclosure trigger and label
##'
##' @param text Label text
##' @param markup logical. Does text have markup? (Toolkit dependent: only implemented for \code{RGtk2}, in \code{qtbase} one can pass HTML formatted text)
##' @param horizontal horizontal (\code{TRUE}) or vertical packing.
##' @param handler handler called when state is toggled
##' @param action passed to handler
##' @param container parent container
##' @param ... passed to parent's \code{add} method
##' @param toolkit toolkit
##' @export
##' @seealso \code{\link{ggroup}} and \code{\link{gframe}}
##' @return An object of class \code{GExpandGroup}
##' inheriting from \code{GFrame}
##' @example inst/examples/ex-gexpandgroup.R

gexpandgroup <- function(
                         text = "", markup = FALSE, horizontal=TRUE,
                         handler = NULL, action = NULL,
                         container = NULL, ... ,
                         toolkit=guiToolkit()){
  obj <- .gexpandgroup (toolkit,
                        text=text, markup=markup, horizontal=horizontal,
                        handler=handler, action=action, container=container ,...
                        )
  visible(obj) <- TRUE               # initial state

  check_return_class(obj, "GExpandGroup")
  obj   
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gexpandgroup
.gexpandgroup <- 
  function(toolkit,
           text = "", markup = FALSE,horizontal=TRUE,
           handler = NULL, action = NULL,
           container = NULL, ... )
  UseMethod( '.gexpandgroup' )


##' visible
##'
##'
##' For gexpandgroup, the visible assignment method is overridden to change the disclosure state
##' @param value logical. If \code{TRUE} show, \code{FALSE} hide.
##' @export
##' @rdname gexpandgroup
##' @usage \method{visible}{GExpandGroup} (obj) <- value
##' @method visible<- GExpandGroup
##' @S3method visible<- GExpandGroup
"visible<-.GExpandGroup" <- function(obj, value) NextMethod()


##' change handler
##'
##' The change handler for a expandGroup is called when the group changes visibility
##' @inheritParams addHandlerChanged
##' @export
##' @rdname gexpandgroup
##' @method addHandlerChanged GExpandGroup
##' @S3method addHandlerChanged GExpandGroup
addHandlerChanged.GExpandGroup <- function(obj, handler, action=NULL, ...) NextMethod()


