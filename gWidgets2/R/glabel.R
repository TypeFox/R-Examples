##' @include methods.R
NULL

##' Basic label widget
##'
##' The basic label widget allows one to label areas of a GUI using
##' text. The most common use would be to label fields in a form. For
##' \pkg{gWidgets2} labels may be editable or responsive to mouse
##' clicks, although it is the author's experience that such uses are
##' not expected by the end user.
##' @param text character. Collapsed using a newline to a single string.
##' @param markup logical. If toolkit supports markup, this indicates
##' it will be used. It is suggested that the \code{font<-} method be
##' used, though for \pkg{gWidgets2Qt} \code{markup} is more
##' convenient.
##' @param editable If TRUE, then clicking on label will enable user-editing of the text.
##' @param handler optional handler. If given, added through addHandlerChanged. Overridden if \code{editable=TRUE}.
##' @param action passed to handler through \code{action} component of first argument of handler. For buttons, this may also be a \code{GAction} instance.
##' @param container parent container (Optional for some toolkits, but not all).
##' @param ... passed to \code{add} method of parent container
##' @param toolkit toolkit instance 
##' @return a \code{GLabel} instance. While this object has its own (reference) methods, one primarily interacts with it through S3 methods defined within the package.
##' @export
##' @examples
##' \dontrun{
##' w <- gwindow("gformlayout", visible=FALSE)
##' g <- gvbox(container=w)
##' g$set_borderwidth(10)
##' 
##' l1 <- glabel("static label", container=g)
##' l2 <- glabel("bold label", container=g)
##' font(l2) <- list(weight="bold")
##' l3 <- glabel("editable label. Click me", editable=TRUE, container=g)
##' 
##' visible(w) <- TRUE
##' 
##' }
glabel <- function(text="", markup=FALSE, editable=FALSE,
                   handler=NULL, action=NULL, container=NULL,
                    ..., toolkit=guiToolkit()) {
  text <- paste(text, collapse="\n")

  if(is.character(toolkit))
    toolkit <- guiToolkit(toolkit)
  
  obj <- .glabel(toolkit, text, markup, editable, handler, action, container, ...)

  check_return_class(obj, "GLabel")
  obj   
  
}

##' S3 generic whose methods are implemented in the toolkit packages
##'
##' @rdname glabel
##' @export
##' @author john verzani
.glabel <- function(toolkit, text, markup=FALSE, editable=FALSE,
                                           handler=NULL, action=NULL, container=NULL,
                                           ...) UseMethod(".glabel")


##' \code{svalue<-} method for a glabel
##'
##' The \code{svalue} methods refer to the main property of the label, its text.
##' @inheritParams svalue<-
##' @export
##' @usage \method{svalue}{GLabel} (obj, index=TRUE, ...) <- value
##' @rdname glabel
##' @method svalue<- GLabel
##' @S3method svalue<- GLabel
"svalue<-.GLabel" <- function(obj, index=TRUE,  ..., value) {
  value <- paste(value, collapse="\n")
  NextMethod()
}


