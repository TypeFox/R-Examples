##' @include methods.R
NULL

##' menu constructor, main and popup
##'
##' A menu may be viewed as a heirarchical collection of buttons, each
##' invoked by clicking on the button. These buttons are exposed
##' through submenus. More generally, a widget may replace the
##' button. This widget intends to support buttons (gactions),
##' separators (gseparator), radio button (gradio) and checkbutton
##' (gcheckbox), but this may be toolkit independent. When using a
##' radio button or checkbox, one should pass in a \code{parent}
##' argument to the constructor -- not a \code{container}.
##' @param menu.list A list defining the menu structure. Named sub
##' lists determine the submenu titles and structure. The list may
##' have components of class: \code{GAction}, mapped to a button;
##' \code{GSeparator}, mapped to a horizontal separator;
##' \code{GRadio}, mapped to linked buttons; or \code{GCheckbox},
##' mapped to a checkbox button.
##' @param popup logical. If true, make a popup window to be added through a handler call
##' @inheritParams gwidget
##' @export
gmenu <- function(
                  menu.list,
                  popup = FALSE,
                  container = NULL,      ... ,
                  toolkit=guiToolkit()){

  deprecated_args <- list(menulist="Use menu.list instead",
                          action="No action argument, parameterize gaction objects individually")
  check_deprecated(deprecated_args, ...)

 
  obj <- .gmenu(toolkit,
                menu.list=menu.list, popup=popup, container=container , ...
                )
  check_return_class(obj, "GMenuBar")
  return(obj)
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gmenu
.gmenu <-  function(toolkit,
                    menu.list=list(),
                    popup=FALSE,
                    container = NULL,
                    ... )
  UseMethod( '.gmenu' )

##' add menubar items to a menu
##'
##' @inheritParams add
##' @param child list. a menubar list or gmenu instance.
##' @export
##' @rdname gmenu
##' @method add GMenuBar
##' @S3method add GMenuBar
add.GMenuBar <- function(obj, child, expand=FALSE, fill=NULL, anchor=NULL, ...) {
  dispatcher <- function(obj, child) UseMethod("dispatcher")
  dispatcher.GMenuBar <- function(child, obj) obj$add_menu_items(obj$widget, svalue(child))
  dispatcher.list <- function(obj, child) obj$add_menu_items(obj$widget, child)
  dispatcher(child, obj)
}


##' "svalue" method
##'
##' For a menubar, \code{svalue} returns the list of action items
##' etc. that defined the menubar. This can be useful to access the
##' underlying item being proxied. (For \code{gaction} items the
##' \code{enabled<-} method may be used on the item, but this may not
##' extend to \code{gradio} and \code{gcheckbox} items)
##' @param index ignored
##' @param drop ignored
##' @export
##' @rdname gmenu
##' @method svalue GMenuBar
##' @S3method svalue GMenuBar
"svalue.GMenuBar" <- function(obj, index=NULL, drop=NULL,...) NextMethod()

##' "svalue<-" method
##'
##' for a menubar, \code{svalue<-} replaces the menubar items with new
##' ones specified by value.
##' @param value a list or menu bar specifying the new men ubar
##' @export
##' @usage \method{svalue}{GMenuBar} (obj, index=NULL, ...) <- value
##' @rdname gmenu
##' @method svalue<- GMenuBar
##' @S3method svalue<- GMenuBar
"svalue<-.GMenuBar" <- function(obj, index=NULL, ..., value) NextMethod()
