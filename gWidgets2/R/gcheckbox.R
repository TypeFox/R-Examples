##' @include methods.R
NULL

##' constructor for checkbox widget
##'
##' A checkbox widget is used to toggle the state of a labeled boolean
##' variable. The main property of this widget is that state, not the
##' label. This variable may be proxied in the usual way -- with a box
##' that indicates or check if \code{TRUE} -- or through a toggle
##' button. 
##' @param text label text
##' @param checked is button selected
##' @param use.togglebutton Use a toggle button (shows depressed) not a check box
##' @param handler Callback called when toggle is changed.
##' @param action passed to handler
##' @param container parent container
##' @param ... passed to \code{add} method of container
##' @param toolkit toolkit
##' @export
##' @return Returns an object of class \code{GCheckbox}.
##' @example inst/examples/ex-selectionwidgets.R
gcheckbox <- function(
                      text="", checked = FALSE, use.togglebutton=FALSE,
                      handler = NULL, action = NULL, container = NULL, ... ,
                      toolkit=guiToolkit()){

  ## text is just first value
  text <- as.character(text)[1]
  ## checked is logical
  checked <- as.logical(checked)[1]
  
  obj <- .gcheckbox (toolkit,
                     text=text, checked=checked,
                     use.togglebutton=use.togglebutton,
                     handler=handler, action=action, container=container, ...
                     )

  check_return_class(obj, "GCheckbox")
  obj
  
}


##' Generic for toolkit dispatch
##'
##' @export
##' @rdname gcheckbox
.gcheckbox <- function(toolkit,
                       text, checked = FALSE, use.togglebutton=FALSE, handler = NULL, action = NULL,
                       container = NULL, ... ) UseMethod( '.gcheckbox' )




##' The change handler for a gcheckbox
##'
##' The change handler for \code{GCheckbox} is called when the value
##' toggles. You can inpsect the current value in the callback to have
##' an action based on the state.
##'
##' @param obj receiver object
##' @export
##' @rdname gcheckbox
##' @method addHandlerChanged GCheckbox
##' @S3method addHandlerChanged GCheckbox
addHandlerChanged.GCheckbox <- function(obj, handler, action=NULL, ...) NextMethod()


##' svalue method
##'
##' The object state is referred to by svalue as a logical (TRUE for checked).
##' The \code{svalue<-} method ensures the value is a logical vector
##' of length 1.
##' @param index ignored. Input is coerced to logical.
##' @param value assignment value
##' @export
##' @usage \method{svalue}{GCheckbox} (obj, index=NULL, ...) <- value
##' @rdname gcheckbox
##' @method svalue<- GCheckbox
##' @S3method svalue<- GCheckbox
"svalue<-.GCheckbox" <- function(obj, index=NULL,  ...,value) {
  value <- as.logical(value)[1]
  NextMethod()
}


##' items assignment takes string
##'
##' The item to select is referred to by the \code{[} method, with only the first element being used.
##' @param x checkbox object
##' @param i item index
##' @param j ignored
##' @note The value is coerced to character, then only first element
##' used for checkbox label
##' @export
##' @usage \method{[}{GCheckbox} (x, i, j, ...) <- value
##' @rdname gcheckbox
##' @method [<- GCheckbox
##' @S3method [<- GCheckbox
"[<-.GCheckbox" <- function(x, i, j, ..., value) {
  value <- as.character(value)[1]
  NextMethod()
}
