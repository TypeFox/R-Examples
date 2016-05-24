##' @include methods.R
NULL


##' A toolbar constructor
##'
##' A toolbar can be added to a main window to proxy various
##' actions. Toolbars can also contain various widgets, such as
##' buttons, checkboxes, radio buttons, etc. These should be
##' constructed using a \code{parent} argument -- not a
##' \code{container} argument. In \pkg{gWidgets2} a toolbar is
##' specified by a list of toolbar items. The \code{svalue} and
##' \code{svalue<-} methods may be used to get or set the items.
##' @param toolbar.list list. A one-level list of \code{gaction}
##' items, \code{gseparator} items or possibly other widgets. In the
##' latter cases the \code{container} argument is not specified
##' prior. (XXX Need to work this out with gWidgetstcltk)
##' @param style style for icon or text.
##' @param container a GWindow instance
##' @param ... ignored
##' @param toolkit toolkit
##' @export
gtoolbar <- function(
                     toolbar.list=list(),
                     style = c("both", "icons", "text", "both-horiz"),
                     container = NULL,
                     ... ,
                     toolkit=guiToolkit()){

  deprecated_args <- list(toolbarlist="Use toolbar.list instead",
                          action="No action argument, parameterize gaction objects individually")
  check_deprecated(deprecated_args, ...)
    
  obj <- .gtoolbar (toolkit,
                    toolbar.list=toolbar.list,
                    style=match.arg(style),
                    container=container ,...
                    )
  check_return_class(obj, "GToolBar")
  return(obj)
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gtoolbar
.gtoolbar <- function(toolkit,
                      toolbar.list=list(),
                      style = c("both", "icons", "text", "both-horiz"),
                      container = NULL,
                      ... )
           UseMethod( '.gtoolbar' )


##' add toolbar items to toolbar
##'
##' A toolbar item is a list of action items or a toolbar instance
##' @inheritParams add
##' @export
##' @rdname gtoolbar
##' @method add GToolBar
##' @S3method add GToolBar
add.GToolBar <- function(obj, child, expand=FALSE, fill=NULL, anchor=NULL, ...) {
  dispatcher <- function(obj, child) UseMethod("dispatcher")
  dispatcher.GToolBar <- function(child, obj) obj$add_toolbar_items(svalue(child))
  dispatcher.list <- function(obj, child) obj$add_toolbar_items(child)
  dispatcher(child, obj)
}




##' "svalue<-" method
##'
##' for a toolbar, \code{svalue<-} replaces the toolbar items with new ones specified by value.
##' @inheritParams svalue
##' @usage \method{svalue}{GToolBar} (obj, index=NULL, ...) <- value
##' @export
##' @rdname gtoolbar
##' @method svalue<- GToolBar
##' @S3method svalue<- GToolBar
"svalue<-.GToolBar" <- function(obj, index=NULL, ..., value) NextMethod()
