##' @include methods.R
NULL

##' Basic button widget
##'
##' The basic button widget is a standard means to provide the user a
##' mechanism to invoke an action. This action may be specified by a
##' handler or by a \code{gaction} object. The main property for
##' \code{GButton} is the label text. If this text matches a stock
##' icon name and the toolkit supports it, an icon will accompany the
##' button.
##' 
##' @param text label text. If text matches a stock icon name, that is used as well
##' @inheritParams gwidget
##' @return a \code{GButton} instance. While this object has its own
##' (reference) methods, one primarily interacts with it through S3
##' methods defined within the package.
##' @export
##' @example inst/examples/ex-gbutton.R
gbutton <- function(text="",   handler=NULL, action=NULL, container=NULL, ..., toolkit=guiToolkit()) {

  if(is.character(toolkit))
    toolkit <- guiToolkit(toolkit)

  deprecated_args <- list(border=c("Border argument has been deprecated.","See reference class method remove_border for availability"))
  check_deprecated(deprecated_args, ...)

  obj <- .gbutton(toolkit, text, handler, action, container, ...)

  check_return_class(obj, "GButton")
  obj
  
}

##' S3 generic whose methods are implemented in the toolkit packages
##'
##' @rdname gbutton
##' @export
.gbutton <- function(toolkit, text, handler, action, container, ...) UseMethod(".gbutton")

##' Change handler for a button is called when a button is activated (e.g., a mouse click)
##'
##' @inheritParams addHandlerChanged
##' @export
##' @rdname gbutton
##' @method addHandlerChanged GButton
##' @S3method addHandlerChanged GButton
addHandlerChanged.GButton <- function(obj, handler, action=NULL, ...) NextMethod()



##' svalue method
##'
##' The \code{svalue} method for a button object refers to its main property, the button label
##' @inheritParams svalue
##' @export
##' @rdname gbutton
##' @method svalue GButton
##' @S3method svalue GButton
svalue.GButton <- function(obj, index=NULL, drop=NULL, ...)   NextMethod()

