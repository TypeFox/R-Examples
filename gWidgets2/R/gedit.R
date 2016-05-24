##' @include methods.R
NULL


##' Single line text edit constructor
##'
##' The default change handler is called when the return key is
##' pressed. It can be useful to also call a handler when the widget
##' loses focus. For that, the \code{addHandlerBlur} method is of
##' use. (This was the default, but is now not, as it was hard to
##' decouple the two when that was desirable.)
##' @param text initial text
##' @param width number of characters
##' @param coerce.with A function or name of function to coerce value with before returning by \code{svalue}
##' @param initial.msg If no initial text is given but an initial
##' message is, then this message is displayed until the widget
##' receives the focus
##' @param handler Change handler. Called when return key is hit. Use
##' \code{addHandleBlur} to add a handler when the widget loses focus,
##' such as through tab-key navigation.
##' @param action passed to handler
##' @param container parent container
##' @param ... passed to \code{add} method of parent
##' @param toolkit toolkit
##' @return An object of class \code{GEdit}. This has sub-classed methods:
##'
##' \enumerate{
##'
##' \item \code{}
##'
##' \item \code{svalue} to retrieve the text
##'
##' \item \code{svalue<-} to set the text
##'
##' \item \code{[} to get the autocomplete values
##'
##' \item \code{[<-} Character. To set autocomplete values
##'
##' \item \code{visible<-} to specify a character to display instead of text (for passwords)
##'
##' }
##'
##' @export
gedit <- function(
                  text = "", width = 25, coerce.with = NULL, initial.msg="",
                  handler = NULL, action = NULL, container = NULL, ... ,
                  toolkit=guiToolkit()) {

  obj <- .gedit(toolkit,
                text=text, width=width, coerce.with=coerce.with, initial.msg=initial.msg,
                handler=handler, action=action, container=container ,...
                )

  check_return_class(obj, "GEdit")
  obj   
}


##' generic for toolkit dispatch
##' 
##' @rdname gedit
##' @export
.gedit <-  function(toolkit,
                    text = "", width = 25, coerce.with = NULL, initial.msg="",
                    handler = NULL, action = NULL, container = NULL, ... ) {
  UseMethod( '.gedit' )
}



##' change handler
##' 
##' The default change handler call is when the user activates the
##' entry by pressing the enter key. Other possible events to consider are
##' covered by: \code{addhandlerBlur} (when the widget loses focuses)
##' and \code{addHandlerKeystroke} (called after each keystroke). For
##' the latter, if the toolkit supports it, the handler's first
##' argument has a component \code{key} passing back the keystroke
##' information.
##' @inheritParams addHandlerChanged
##' @export
##' @rdname gedit
##' @method addHandlerChanged GEdit
##' @S3method addHandlerChanged GEdit
addHandlerChanged.GEdit <- function(obj, handler, action=NULL, ...) NextMethod()



##' svalue method
##'
##' The \code{svalue} method for a edit object refers to its main property, the text in the box. 
##' @inheritParams svalue
##' @export
##' @rdname gedit
##' @method svalue GEdit
##' @S3method svalue GEdit
svalue.GEdit <- function(obj, index=NULL, drop=NULL, ...)   NextMethod()


##' Set available words for autocompletion
##'
##' The underlying widget may allow autocompletion, if this is the
##' case then this method is used to set the list of candidates.
##' @inheritParams gedit
##' @export
##' @rdname gWidgets2-S3methods
##' @method [ GEdit
##' @S3method [ GEdit
"[.GEdit" <- function(x, i, j, ..., drop=TRUE) NextMethod()

