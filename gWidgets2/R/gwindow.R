##' @include methods.R
NULL


##' top-level window object
##'
##' @title gwindow
##' @param title title for window's title bar. This is the main
##' property and is accessed via \code{svalue} or \code{svalue<-}.
##' @param visible logical. If code{TRUE} window is drawn when
##' constructed. Otherwise, window can be drawn later using
##' \code{visible<-}. This value can default to \code{FALSE} by
##' setting the option:
##' \code{options("gWidgets:gwindow-default-visible-is-false"=TRUE)}.
##' There are advantages: windows can draw slowly when adding many
##' items. With \pkg{gWidgets2RGtk2}, the  \code{ggraphics} widget can
##' like to be added to an undrawn widget as this avoids sizing issue.
##' @param name Name for registry of windows
##' @param width initial width of window
##' @param height initial height of window. This sets height before window manager manages the window
##' @param parent If non-NULL, can be used to suggest default location
##' of window. The argument name was changed from location to
##' parent. This can be a coordinate pair (x,y) with (0,0) the upper
##' left corner, or a gwindow instance. In the latter case the
##' location is suggested by the location of the current window. This
##' is useful for placing dialogs near the parent window.
##' @param handler handler for destroy event
##' @param action action passed t handler
##' @param ... ignored
##' @param toolkit toolkit
##' @return a \code{GWindow} instance
##' @export
gwindow <- function(title="Window", visible=TRUE, name=title, width=NULL, height=NULL, parent=NULL, handler=NULL, action=NULL, ..., toolkit=guiToolkit()) {

  ## process arguments
  deprecated_args=list(
    "location", "This argument was renamed 'parent'")
  check_deprecated(...)

  theArgs <- list(...)
  if(!is.null(theArgs$location)) {
    parent <- theArgs$location
  }
  
  ## THe visible=TRUE default is not the best. I'd change it if I could go back in time, but
  ## c'est la vie. Anyways, for those that it really bugs there is this check
  if(!is.null(getOption("gWidgets:gwindow-default-visible-is-false")))
    visible <- FALSE

  if(is.character(toolkit))
    toolkit <- guiToolkit(toolkit)
  
  ## now go forth and dispatch
  obj <- .gwindow(toolkit, title, visible=visible, name, width, height, parent, handler, action,  ...)
  check_return_class(obj, "GWindow")
  obj
}

##' S3 generic whose methods are implemented in the toolkit packages
##'
##' @rdname gwindow
##' @export
##' @author john verzani
.gwindow <- function(toolkit, title, visible, name, width, height, parent, handler , action, ...) UseMethod(".gwindow")



##' add method for top-level windows
##'
##' Dispatches on type of child (menubar, toolbar, statusbar, widget)
##' @inheritParams add
##' @export
##' @rdname gwindow
##' @method add GWindow
##' @S3method add GWindow
add.GWindow <- function(obj, child, expand=NULL, fill=NULL, anchor=NULL, ...) {
  ## one way -- poor man's way -- of double dispatch on obj and class of child
  .add <- function(child, obj, ...) UseMethod(".add")
  .add.default <- function(child, obj, ...) obj$add_child(child)
  .add.GMenu <- function(child, obj, ...) obj$add_menu(child)
  .add.GToolBar <- function(child, obj, ...) obj$add_toolbar(child)
  .add.GStatusbar <- function(child, obj, ...) obj$add_statusbar(child)
  
  .add(child, obj, ...)
}
  
##' dispose method for gwindow
##'
##' The \code{dispose} method destroys the top-level window and its children.
##' @inheritParams dispose
##' @export
##' @rdname gwindow
##' @method dispose GWindow
##' @S3method dispose GWindow
dispose.GWindow <- function(obj, ...) {
  obj$dispose_window()
}
  
