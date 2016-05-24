##' @include BasicInterface.R
NULL


## Handlers


##' Generic method to add a handler passing a signal
##'
##' A GUI is made interactive by assigning handlers to user-generated
##' events, such as a mouse click, change of widget state, or keyboard
##' press. In \pkg{gWidgets2} handlers are assigned through various
##' \code{addHandlerXXX} methods. The handlers are functions whose
##' first argument should expect a list with components \code{obj} (to
##' pass in the receiver object) and \code{action} (to pass in any
##' user-supplied value to the \code{action} argument). Some handlers
##' add other components, such as mouse position information on a
##' click, or key information on a keyboard event.
##'
##' Although the \code{add_handler} method, to which \code{addHandler}
##' dispatches, is basically the workhorse to add a handler to
##' response to a signal, it generally isn't called directly, as its
##' use is not cross toolkit. Rather, if possible, one should use the
##' \code{addHandlerXXX} methods to add a handler. These dispatch to
##' this (basically) but do so in a toolkit independent manner.
##'
##' This call (and the others) returns a handler ID which may be used
##' for some toolkits later on to remove, block or unblock the
##' call. All handlers for a widget may be blocked or unblocked via
##' \code{blockHandlers} and \code{unblockHandlers}.
##' 
##' @param obj object receiving event and emitting a signal to the handler
##' @param signal toolkit signal, e.g. "clicked"
##' @param handler handler to assign when signal is emitted. A handler
##' is a function, its first argument should expect a list with
##' components \code{obj} containing a reference to the object and
##' \code{action}. Some handlers are passed additional values.
##' @param action passed to handler to parameterize call.
##' @param ... passed along
##' @note This method is not toolkit independent, as the signal value depends on the toolkit
##' @return a handler ID which can be used to block/unblock or remove the handler
##' @seealso \code{\link{blockHandlers}},
##' \code{\link{unblockHandlers}}, \code{\link{blockHandler}},
##' \code{\link{unblockHandler}}, and \code{\link{removeHandler}}
##' @export
##' @rdname gWidgets-handlers
addHandler <- function(obj, signal, handler, action=NULL, ...) UseMethod("addHandler")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandler default
##' @S3method addHandler default
addHandler.default <- function(obj, signal, handler, action=NULL, ...) 
  obj$add_handler(signal, handler, action=action, ...) 

##' Add a handler to the generic "changed" event, which is the main event for a widget
##'
##' The "changed" event varies wildly amongst the widgets, but is
##' meant to be the most "obvious" one. Typically this is also similar
##' to "selected".
##'
##' The "changed" event is also the one that a handler
##' passed to the constructor is called on. 
##' @export
##' @rdname gWidgets-handlers
addHandlerChanged <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerChanged")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerChanged default
##' @S3method addHandlerChanged default
addHandlerChanged.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_changed(handler, action=action, ...)

##' Add handler for clicked event
##'
##' @export
##' @rdname gWidgets-handlers
addHandlerClicked <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerClicked")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerClicked default
##' @S3method addHandlerClicked default
addHandlerClicked.default <-  function(obj, handler, action=NULL, ...)
    obj$add_handler_clicked(handler, action=action, ...)

##' Add handler for double click event
##'
##' @export
##' @rdname gWidgets-handlers
addHandlerDoubleclick <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerDoubleclick")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerDoubleclick default
##' @S3method addHandlerDoubleclick default
addHandlerDoubleclick.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_double_clicked(handler, action=action, ...)


##' Add handler for right click event
##'
##' This may not be supported by all toolkits.
##' @export
##' @rdname gWidgets-handlers
addHandlerRightclick <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerRightclick")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerRightclick default
##' @S3method addHandlerRightclick default
addHandlerRightclick.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_right_clicked(handler, action=action, ...)


##' Add handler for shift click event
##'
##' This may not be supported by all toolkits.
##' @export
##' @rdname gWidgets-handlers
addHandlerShiftclick <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerShiftclick")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerShiftclick default
##' @S3method addHandlerShiftclick default
addHandlerShiftclick.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_shift_clicked(handler, action=action, ...)


##' Add handler for control+click event
##'
##' This may not be supported by all toolkits.
##' @export
##' @rdname gWidgets-handlers
addHandlerControlclick <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerControlclick")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerControlclick default
##' @S3method addHandlerControlclick default
addHandlerControlclick.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_control_clicked(handler, action=action, ...)




##' Add handler for column click event
##'
##' For table widgets (\code{gtable}, \code{gdf}) clicking the column
##' header should trigger this event. The column that is clicked on is
##' passed back in the component \code{column}.
##' @export
##' @rdname gWidgets-handlers
addHandlerColumnclicked <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerColumnclicked")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerColumnclicked default
##' @S3method addHandlerColumnclicked default
addHandlerColumnclicked.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_column_clicked(handler, action=action, ...)

##' Add handler for column double click event
##'
##' If defined (\code{gtable}, \code{gdf}) calls event handler for
##' double click enent. Passes back column information through
##' \code{column} component.
##' @export
##' @rdname gWidgets-handlers
addHandlerColumnDoubleclicked <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerColumnDoubleclicked")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerColumnDoubleclicked default
##' @S3method addHandlerColumnDoubleclicked default
addHandlerColumnDoubleclicked.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_column_double_clicked(handler, action=action, ...)


##' Add handler for column right click event
##'
##' @export
##' @rdname gWidgets-handlers
addHandlerColumnRightclicked <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerColumnRightclicked")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerColumnRightclicked default
##' @S3method addHandlerColumnRightclicked default
addHandlerColumnRightclicked.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_column_right_clicked(handler, action=action, ...)



##' Add a handler to the a "select" event
##'
##' The select event defaults to the "changed" event.
##' @export
##' @rdname gWidgets-handlers
addHandlerSelect <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerSelect")

##' Default S3 method
##'
##' @inheritParams addHandlerSelect
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerSelect default
##' @S3method addHandlerSelect default
addHandlerSelect.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_select(handler, action=action, ...)



##' Add a handler to the a "selection-changed" event
##'
##' The "select" event is when a user "selects" an object, the
##' "selection changed" event is when the selection changes. The
##' distinction is in table and tree widgets where a user may select
##' values with a single click yet want to initiate an action with a
##' double click. The latter is the "addHandlerSelect" event, the
##' former the "addHandlerSelectionChanged" event.
##' @export
##' @rdname gWidgets-handlers
addHandlerSelectionChanged <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerSelectionChanged")

##' Default S3 method
##'
##' @inheritParams addHandlerSelectionChanged
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerSelectionChanged default
##' @S3method addHandlerSelectionChanged default
addHandlerSelectionChanged.default <- function(obj, handler, action=NULL, ...)
  obj$add_handler_selection_changed(handler, action=action, ...)


##' Add handler for focus in event
##'
##' When a widget has the focus, it will receive the keyboard
##' input. This handler is called when a widget receives the focus.
##' @export
##' @rdname gWidgets-handlers
addHandlerFocus <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerFocus")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerFocus default
##' @S3method addHandlerFocus default
addHandlerFocus.default <-   function(obj, handler, action=NULL, ...)
  obj$add_handler_focus(handler, action=action, ...)

##' Add handler for blur, or focus-out, event
##'
##' A blur or focus out event for a widget triggers this event handler
##' @export
##' @rdname gWidgets-handlers
addHandlerBlur <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerBlur")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerBlur default
##' @S3method addHandlerBlur default
addHandlerBlur.default <-  function(obj, handler, action=NULL, ...)
  obj$add_handler_blur(handler, action=action, ...)


##' Add handler for destroy event
##'
##' When a widget is destroyed, a handler can be assigned to perform any clean up tasks that are needed.
##' @seealso \code{\link{addHandlerUnrealize}}.
##' @export
##' @rdname gWidgets-handlers
addHandlerDestroy <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerDestroy")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerDestroy default
##' @S3method addHandlerDestroy default
addHandlerDestroy.default <-  function(obj, handler, action=NULL, ...)
  obj$add_handler_destroy(handler, action=action, ...)

##' Add handler for unrealize event for a top-level window
##'
##' For gwindow objects this handler is called before the window is
##' closed. If this handler returns \code{TRUE} the window will be
##' closed, if \code{FALSE} the window will not be closed. In
##' contrast, the "destroy" handler does not allow conditional
##' destruction.
##' @export
##' @rdname gWidgets-handlers
addHandlerUnrealize <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerUnrealize")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerUnrealize default
##' @S3method addHandlerUnrealize default
addHandlerUnrealize.default <-  function(obj, handler, action=NULL, ...)
  obj$add_handler_unrealize(handler, action=action, ...)





##' Add handler for expose event (when a widget is exposed, say it had been covered)
##'
##' @export
##' @rdname gWidgets-handlers
addHandlerExpose <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerExpose")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerExpose default
##' @S3method addHandlerExpose default
addHandlerExpose.default <-  function(obj, handler, action=NULL, ...)
  obj$add_handler_expose(handler, action=action, ...)



##' Add handler for keystroke events
##'
##' The "h" argument has components \code{key} for the key and possibly \code{modifier} for the modifier.
##' @export
##' @rdname gWidgets-handlers
addHandlerKeystroke <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerKeystroke")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerKeystroke default
##' @S3method addHandlerKeystroke default
addHandlerKeystroke.default <-  function(obj, handler, action=NULL, ...)
  obj$add_handler_keystroke(handler, action=action, ...)



##' Add handler for mousemotion events
##'
##' @export
##' @rdname gWidgets-handlers
addHandlerMouseMotion <- function(obj, handler, action=NULL, ...) UseMethod("addHandlerMouseMotion")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addHandlerMouseMotion default
##' @S3method addHandlerMouseMotion default
addHandlerMouseMotion.default <-  function(obj, handler, action=NULL, ...)
  obj$add_handler_mouse_motion(handler, action=action, ...)

##' Add an idle handler
##'
##' deprecated. See \code{\link{gtimer}}.
##' @export
##' @rdname gWidgets-handlers
addHandlerIdle <- function( ...) {
  .Deprecated("No addHandleIdle method. Use gtimer for that purpose")
}

##' Add a "popup" menu to the widget
##'
##' Defaults to adding a right-click mouse popup menu, better known as a
##' context menu, though some toolkits have both this and the latter
##' provided. 
##' @param menulist a list of \code{gaction} items or a \code{gmenu} instance
##' @export
##' @rdname gWidgets-handlers
addPopupMenu <- function(obj, menulist, action=NULL, ...) UseMethod("addPopupMenu")

##' S3 method for popup menu
##'
##' @export
##' @rdname gWidgets-handlers
##' @method addPopupMenu default
##' @S3method addPopupMenu default
addPopupMenu.default <-  function(obj, menulist, action=NULL, ...)
  obj$add_popup_menu(menulist, action=action, ...)


##' Add a context  "popup" menu to the widget
##'
##' These menus are also known as context menus, though there isn't
##' really a good mechanism within \pkg{gWidgets2} to make the menu
##' items context sensitive.
##' @inheritParams addPopupMenu
##' @export
##' @rdname gWidgets-handlers
addRightclickPopupMenu <- function(obj, menulist, action=NULL, ...) UseMethod("addRightclickPopupMenu")

##' S3 method for popup menu
##'
##' @export
##' @rdname gWidgets-handlers
##' @method addRightclickPopupMenu default
##' @S3method addRightclickPopupMenu default
addRightclickPopupMenu.default <-  function(obj, menulist, action=NULL, ...)
  obj$add_3rd_mouse_popup_menu(menulist, action=action, ...)


##' Alias for poor initial name choice
##'
##' @rdname gWidgets-handlers
##' @export
##' @aliases addRightclickPopupMenu
add3rdmousePopupMenu <- addRightclickPopupMenu

##' Alias for poor punctation choice
##'
##' @rdname gWidgets-handlers
##' @export
##' @aliases addRightclickPopupMenu
add3rdMousePopupMenu <- add3rdmousePopupMenu

##' Specify a widget is a source for a drop action
##'
##' Drag and drop requires one to register widgets a sources for
##' dragging, a widgets as a targets for dropping.
##'
##' To specify the values that is transferred in a drag and drop
##' event, the handler specified here should return the value to pass
##' via drag and drop. It will appear as the \code{dropdata} component
##' of the list passed in as the first argument of the drop handler
##' @inheritParams addHandler
##' @param data.type Type of data returned. It is either text or an object
##' @export
##' @rdname gWidgets-handlers
addDropSource <- function(obj, handler, action=NULL, data.type=c("text", "object"), ...) UseMethod("addDropSource")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addDropSource default
##' @S3method addDropSource default
addDropSource.default <-  function(obj, handler, action=NULL, data.type=c("text", "object"), ...)
  obj$add_drop_source(handler, action=action, data.type=match.arg(data.type), ...)

##' Specify that a widget is a drop target
##'
##' The handler is called on the drop event. The component
##' \code{dropdata} passes in the value being transferred by dragging.
##' @export
##' @rdname gWidgets-handlers
addDropTarget <- function(obj, handler, action=NULL, ...) UseMethod("addDropTarget")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addDropTarget default
##' @S3method addDropTarget default
addDropTarget.default <-  function(obj, handler, action=NULL, ...)
  obj$add_drop_target( handler, action=action, ...)


##' When a drag event crosses over the object the handler is called. 
##'
##' @export
##' @rdname gWidgets-handlers
addDragMotion <- function(obj, handler, action=NULL, ...) UseMethod("addDragMotion")

##' Default S3 method
##'
##' @inheritParams addHandler
##' @export
##' @rdname gWidgets-handlers
##' @method addDragMotion default
##' @S3method addDragMotion default
addDragMotion.default <-  function(obj, handler, action=NULL, ...)
  obj$add_drag_motion( handler, action=action, ...)

##' block all handlers for object
##'
##' Block all handlers for an object. Removed via unblockHandlers.
##' @return NULL
##' @export
##' @rdname gWidgets-handlers
blockHandlers <- function(obj, ...) UseMethod("blockHandlers")



##' S3 method to block all handlers
##'
##' @export
##' @rdname gWidgets-handlers
##' @method blockHandlers default
##' @S3method blockHandlers default
blockHandlers.default <- function(obj, ...) obj$block_handlers(...)


##' Block a handler
##'
##' @param ID returned by addHandler. If missing will try to block all handler passed to constructor
##' @note For the gWidgets2Qt package one can not block, unblock or
##' remove a single handler, but rather must do all the objects
##' handlers at once.
##' @seealso \code{\link{blockHandlers}} to block all handlers for widget
##' @export
##' @rdname gWidgets-handlers
blockHandler <- function(obj, ID, ...) UseMethod("blockHandler")

##' S3 method to block handler
##'
##' @export
##' @rdname gWidgets-handlers
##' @method blockHandler default
##' @S3method blockHandler default
blockHandler.default <- function(obj, ID, ...) obj$block_handler(ID)



##' method call to unblock global handler block.
##'
##' The block is a counter that gets decremented. If more
##' blockHandlers calls are made than unblockHandlers, the handlers
##' will still be blocked.
##' @export
##' @rdname gWidgets-handlers
unblockHandlers <- function(obj, ...) UseMethod("unblockHandlers")


##' S3 method to block handler
##'
##' @export
##' @rdname gWidgets-handlers
##' @method unblockHandlers default
##' @S3method unblockHandlers default
unblockHandlers.default <- function(obj, ...) obj$unblock_handlers()

##' method call to unblock a blocked handler
##'
##' @export
##' @rdname gWidgets-handlers
unblockHandler <- function(obj, ID, ...) UseMethod("unblockHandler")


##' S3 method to block handler
##'
##' @export
##' @rdname gWidgets-handlers
##' @method unblockHandler default
##' @S3method unblockHandler default
unblockHandler.default <- function(obj, ID, ...) obj$unblock_handler(ID)


##' method call to unblock a remove permanently handler
##'
##' @export
##' @rdname gWidgets-handlers
removeHandler <- function(obj, ID, ...) UseMethod("removeHandler")

##' S3 method to remove handler
##'
##' @export
##' @rdname gWidgets-handlers
##' @method removeHandler default
##' @S3method removeHandler default
removeHandler.default <- function(obj, ID, ...) obj$remove_handler(ID)
