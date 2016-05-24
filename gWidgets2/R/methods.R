##' @include S3-methods.R
NULL

## generic methods and definitions

##' Common parts of a widget
##'
##' Used as template for documentation
##' @param handler A handler assigned to the default change
##' signal. Handlers are called when some event triggers a widget to
##' emit a signal. For each widget some default signal is assumed, and
##' handlers may be assigned to that through \code{addHandlerChanged}
##' or at construction time. Handlers are functions whose first
##' argument, \code{h} in the documentation, is a list with atleast
##' two components \code{obj}, referring to the object emitting the
##' signal and \code{action}, which passes in user-specified data to
##' parameterize the function call.
##'
##' Handlers may also be added via \code{addHandlerXXX} methods for
##' the widgets, where \code{XXX} indicates the signal, with a default
##' signal mapped to \code{addHandlerChanged}
##' (cf. \code{\link{addHandler}} for a listing). These methods pass
##' back a handler ID that can be used with \code{blockHandler} and
##' \code{unblockHandler} to suppress temporarily the calling of the
##' handler.
##' @param action User supplied data passed to the handler when it is called
##' @param container A parent container. When a widget is created it can be
##' incorporated into the widget heirarchy by passing in a parent
##' container at construction time. (For some toolkits this is not
##' optional, e.g. \pkg{gWidgets2tcltk} or \pkg{gWidgets2WWW2}.)
##' @param ... These values are passed to the \code{add} method of the
##' parent container. Examples of values are \code{expand},
##' \code{fill}, and \code{anchor}, although they're not always
##' supported by a given widget. For more details see \link{add}.
##' Occasionally the variable arguments feature has been used to sneak
##' in hidden arguments to toolkit implementations. For example, when
##' using a widget as a menubar object one can specify a parent
##' argument to pass in parent information, similar to how the
##' argument is used with gaction and the dialogs.
##' @param toolkit Each widget constructor is passed in the toolkit it
##' will use. This is typically done using the default, which will
##' lookup the toolkit through \code{\link{guiToolkit}}.
gwidget <- function(handler=NULL, action=NULL, container=NULL, ...,toolkit=guiToolkit()) {}


##' Common parts of a container widget
##'
##' Used as template for documentation
##' @param container A parent container. When a widget is created it can be
##' incorporated into the widget heirarchy by passing in a parent
##' container at construction time. (For some toolkits this is not
##' optional, e.g. \pkg{gWidgets2tcltk} or \pkg{gWidgets2WWW2}.)
##' @param ... These values are passed to the \code{add} method of the
##' parent container, and occasionally have been used to sneak in
##' hidden arguments to toolkit implementations.
##' @param toolkit Each widget constructor is passed in the toolkit it
##' will use. This is typically done using the default, which will
##' lookup the toolkit through \code{\link{guiToolkit}}.
gcontainer <- function(container=NULL, ...,toolkit=guiToolkit()) {}

##' svalue
##'
##' This returns the "selected" value in a widget (if applicable), or
##' its main property. Selection varies from widget to widget, but
##' should generally be what can be added to the widget by mouse click
##' or typing. For some widgets, the extra argument \code{index=TRUE}
##' will return the index of the selected value, not the value. For
##' some widget, the argument \code{drop} is given to either prevent
##' or encourage dropping of information.
##' 
##' @param obj object of method call
##' @param index NULL or logical. If \code{TRUE} and widget supports it an index, instead of a value will be returned.
##' @param drop NULL or logical. If widget supports it, drop will work as it does in a data frame or perhaps someother means.
##' @param ... passed on to call
##' @return THe return value varies, depending if the widget is a
##' "selection" widget or not. For non-selection widgets, the main
##' property is loosely defined (the title of a window, text of a
##' label or button, spacing of box containers, ...). For selection
##' widgets the return value is the currently selected value. If no
##' selection is made, this will be a 0-length vector with the
##' expected class, if possible. For selection widgets, when
##' \code{index=TRUE}, the value is an integer, possible 0-length when
##' non selection is made.
##' @rdname svalue
##' @export
svalue <- function(obj, index=FALSE, drop=NULL, ...) UseMethod("svalue")

##' default svalue instance
##'
##' Calls \code{coerce_with} when available. This value is a function
##' and may be set as any property if the constructor does not
##' explicity provide it.
##' @export
##' @rdname svalue
##' @S3method svalue default
##' @method svalue default
svalue.default <- function(obj, index=NULL, drop=NULL, ...) {
  if(!isExtant(obj)) {
    return()
  }
  if(getWithDefault(index, FALSE)) {
    val <- obj$get_index(drop=drop, ...)
  } else {
    val <- obj$get_value(drop=drop, ...)
    if(exists("coerce_with", obj) &&
       !is(obj$coerce_with, "uninitializedField") &&
       !is.null(obj$coerce_with)) {
      if(is.character(obj$coerce_with))
        obj$coerce_with <- get(obj$coerce_with, inherits=TRUE)
      val <- obj$coerce_with(val)
    }
  }
  val
}

##' svalue<-
##'
##' This method sets the selected value of, or main property of the widget.
##' @param value value to assign for selection or property
##' @rdname svalue
##' @export
##' @usage svalue(obj, index=NULL, ...) <- value
"svalue<-" <- function(obj, index=NULL, ..., value) UseMethod("svalue<-")

##' Base S3 method
##'
##' @rdname svalue
##' @export
##' @usage svalue(obj, index=NULL, ...) <- value
##' @S3method svalue<- default
##' @method svalue<- default
"svalue<-.default" <- function(obj, index=NULL, ..., value) {
  if(!isExtant(obj)) {
    return(obj)
  }

  if(getWithDefault(index, FALSE))
    obj$set_index(value, ...)
  else
    obj$set_value(value, ...)
  obj
}

##' enabled
##'
##' A widget is enabled if it is sensitive to user input
##' @param obj object
##' @export
##' @return logical indicating if widget is enabled
##' @rdname enabled
enabled <- function(obj) UseMethod("enabled")

##' base S3 method for enabled.
##'
##' @export
##' @rdname enabled
##' @S3method enabled default
##' @method enabled default
enabled.default <- function(obj) {
  if(isExtant(obj))
    obj$get_enabled()
}
##' Set whether widget is enabled or not
##'
##' @param value logical
##' @return if \code{value} is logical and \code{FALSE} widget will be insensitive to user input and rendered in a muted state.
##' @export
##' @usage enabled(obj) <- value
##' @rdname enabled
"enabled<-" <- function(obj, value) UseMethod("enabled<-")

##' S3 method for setting enabled property
##'
##' @export
##' @usage enabled(obj) <- value
##' @rdname enabled
##' @S3method enabled<- default
##' @method enabled<- default
"enabled<-.default" <- function(obj, value) {
  if(isExtant(obj))
    obj$set_enabled(as.logical(value))
  obj
}

##' Controls whether widget is visible or not
##'
##' For most -- but not all -- widgets, a widget is visible if it is
##' shown. For others, parts of the widget may be controlled by
##' visible. If the former state is desired, simply place widget into
##' a box container.
##' @param obj object
##' @param ... ignored
##' @export
##' @rdname visible
visible <- function(obj, ...) UseMethod("visible")

##' Basic S3 method
##'
##' @export
##' @rdname visible
##' @S3method visible default
##' @method visible default
visible.default <- function(obj, ...) {
  if(isExtant(obj))
    obj$get_visible()
}

##' Set visibility of an object
##'
##' @param value logical. Set visible state.
##' @export
##' @usage visible(obj) <- value
##' @rdname visible
"visible<-" <- function(obj, value) UseMethod("visible<-")

##' Basic S3 method for visible
##'
##' @export
##' @usage visible(obj) <- value
##' @rdname visible
##' @method visible<- default
##' @S3method visible<- default
"visible<-.default" <- function(obj, value) {
  if(isExtant(obj))
    obj$set_visible(as.logical(value))
  obj
}

##' Does widget have focus
##'
##' a widget has focus if it will receive input events
##' @param obj object
##' @export
##' @rdname focus
focus <- function(obj) UseMethod("focus")

##' Basic S3 method
##'
##' @export
##' @rdname focus
##' @method focus default
##' @S3method focus default
focus.default <- function(obj) {
  if(isExtant(obj))
    obj$get_focus()
}

##' Set focus onto object. 
##'
##' For some widgets, this sets user focus (e.g. gedit gets focus for
##' typing). For others, settig the focus calls the raise
##' methods. (for gwindow, it will raise the window)
##' @param value logical. Set focus state.
##' @export
##' @usage focus(obj) <- value
##' @rdname focus
"focus<-" <- function(obj, value) UseMethod("focus<-")

##' Basic S3 method for focus
##'
##' @export
##' @usage focus(obj) <- value
##' @rdname focus
##' @method focus<- default
##' @S3method focus<- default
"focus<-.default" <- function(obj, value) {
  if(isExtant(obj))
    obj$set_focus(as.logical(value))
  obj
}


##' Controls whether widget is editable or not
##'
##' Some widgets may be editable. If possible, the setter method can
##' be used to toggle the state. This method indicates the state.
##' @param obj object
##' @param i index to apply to, when applicable
##' @export
##' @rdname editable
editable <- function(obj, i) UseMethod("editable")

##' Basic S3 method
##'
##' @export
##' @rdname editable
##' @method editable default
##' @S3method editable default
editable.default <- function(obj, i) {
  if(isExtant(obj))
    obj$get_editable(i)
}
##' Set whether an object can be edited
##'
##' @param value logical. Set editable state.
##' @export
##' @usage editable(obj, i) <- value
##' @rdname editable
"editable<-" <- function(obj, i, value) UseMethod("editable<-")

##' Basic S3 method for editable
##'
##' @export
##' @usage editable(obj, i) <- value
##' @rdname editable
##' @method editable<- default
##' @S3method editable<- default
"editable<-.default" <- function(obj, i,  value) {
  if(isExtant(obj))
    obj$set_editable(as.logical(value), i)
  obj
}

##' Returns font specification for widget, if available
##'
##' @param obj object
##' @export
##' @rdname font
font <- function(obj) UseMethod("font")

##' Basic S3 method for font
##'
##' @export
##' @rdname font
##' @method font default
##' @S3method font default
font.default <- function(obj) {
  if(isExtant(obj))
    obj$get_font()
}
##' Set font for a widget
##'
##' @param value The font specification is given in terms of a named
##' vector or list where the names indicate a font attribute and the
##' value a reasonable choice:
##' 
##' \describe{
##' 
##' \item{weight}{c("light", "normal", "medium", "bold", "heavy")}
##'
##' \item{style}{c("normal", "oblique", "italic")}
##'
##' \item{family}{c("sans", "helvetica", "times", "monospace")}
##'
##' \item{size}{an integer, say c(6,8,10,11,12,14,16,18,20, 24,36,72)}
##'
##' \item{color (or foreground)}{One of colors()}
##'
##' \item{background}{One of colors()}
##'
##' \item{scale}{c("xx-large", "x-large",  "large" ,   "medium",   "small",    "x-small",  "xx-small")}
##'
##' }
##'
##' These are from Gtk's font specs, which though fairly standard, may
##' not be totally supported in the other toolkits.
##' @export
##' @usage font(obj) <- value
##' @rdname font
"font<-" <- function(obj, value) UseMethod("font<-")

##' Basic S3 method for setting font
##'
##' @export
##' @usage font(obj) <- value
##' @rdname font
##' @method font<- default
##' @S3method font<- default
"font<-.default" <- function(obj, value) {
  if(isExtant(obj))
    obj$set_font(value)
  obj
}

##' get a persistent attribute for an object
##'
##' @param obj object
##' @param key character. Values are stored by key. If missing, all keys are returned.
##' @export
##' @rdname tag
tag <- function(obj, key) UseMethod("tag")

##' Basic S3 method
##'
##' @export
##' @rdname tag
##' @method tag default
##' @S3method tag default
tag.default <- function(obj, key) {
  if(isExtant(obj))
    obj$get_attr(key)
}

##' set a persistent attribute for an object
##'
##' Unlike \code{attr<-}, this method (essentially) stores the
##' attribute in a reference to the object, not a copy. As such it can
##' be used within function call (handlers) to assign values outside
##' the scope of the function call.
##' @param value to assign to key
##' @export
##' @usage tag(obj, key) <- value
##' @rdname tag
"tag<-" <- function(obj, key, value) UseMethod("tag<-")

##' Basic S3 method
##'
##' @export
##' @usage tag(obj, key) <- value
##' @rdname tag
##' @method tag<- default
##' @S3method tag<- default
"tag<-.default" <- function(obj, key, value) {
  if(isExtant(obj))
    obj$set_attr(key, value)
  obj
}

## XXX add others  size, size<-, 

##' Return size (width and height) of widget
##'
##' @param obj object
##' @rdname size
##' @export
size <- function(obj) UseMethod("size")

##' S3 method for size
##'
##' @export
##' @rdname size
##' @method size default
##' @S3method size default
size.default <- function(obj) {
  if(isExtant(obj))
    obj$get_size()
}

##' Set size of object (width, height)
##'
##' The size is specified in pixels (integers). Some toolkits allow -1 as a default, but not all.
##' @param value size in pixels
##' @export
##' @usage size(obj) <- value
##' @rdname size
"size<-" <- function(obj, value) UseMethod("size<-")

##' S3 method for size
##'
##' @export
##' @usage size(obj) <- value
##' @rdname size
##' @method size<- default
##' @S3method size<- default
"size<-.default" <- function(obj, value) {
  if(isExtant(obj))
    obj$set_size(value)
  obj
}

##' Get a tooltip for the widget
##'
##' @param obj object
##' @export
##' @rdname tooltip
"tooltip" <- function(obj) UseMethod("tooltip")

##' Basic S3 method for tooltip<-
##'
##' @export
##' @rdname tooltip
##' @method tooltip default
##' @S3method tooltip default
"tooltip.default" <- function(obj) {
  if(isExtant(obj))
    obj$get_tooltip()
}

##' Set a tooltip for the widget
##'
##' @param value character tooltip value
##' @export
##' @usage tooltip(obj) <- value
##' @rdname tooltip
"tooltip<-" <- function(obj, value) UseMethod("tooltip<-")

##' Basic S3 method for tooltip<-
##'
##' @export
##' @usage tooltip(obj) <- value
##' @rdname tooltip
##' @method tooltip<- default
##' @S3method tooltip<- default
"tooltip<-.default" <- function(obj, value) {
  if(isExtant(obj))
    obj$set_tooltip(paste(value, collapse="\n"))
  obj
}


##' Undo past action. 
##'
##' Some widgets support undo actions. See reference class method \code{can_undo} as well.
##' @param obj object to call undo on
##' @param ... ignored
##' @export
##' @rdname undo
undo <- function(obj, ...) UseMethod("undo")

##' S3 method. 
##'
##' @export
##' @rdname undo
##' @S3method undo GComponent
##' @method undo GComponent
undo.GComponent <- function(obj, ...) {
  if(isExtant(obj))
    obj$undo(...)
}



##' Redo past action. 
##'
##' Some widgets support redo actions
##' @param obj object to redo
##' @param ... ignored
##' @export
##' @rdname redo
redo <- function(obj, ...) UseMethod("redo")

##' S3 method. 
##'
##' @export
##' @rdname redo
##' @method redo GComponent
##' @S3method redo GComponent
redo.GComponent <- function(obj, ...) {
  if(isExtant(obj))
    obj$redo(...)
}

##' Check if widget is extant.
##'
##' Widgets can be destroyed, but their R object is still present. This is FALSE in that case.
##' @param obj object
##' @export
##' @rdname isExtant
"isExtant" <- function(obj) UseMethod("isExtant")

##' Basic S3 method for isExtant
##'
##' @export
##' @rdname isExtant
##' @method isExtant default
##' @S3method isExtant default
"isExtant.default" <- function(obj) {
  ret <- try(obj$is_extant(), silent=TRUE)
  if(is(ret, "try-error"))
    FALSE
  else
    ret
}


## container methods

##' Add a child object to parent container
##'
##' Add packs in child objects.
##' @param obj parent object
##' @param child child widget
##' @param expand NULL or logical. For box containers controls whether a child will expand to fill the allocated space. 
##' @param fill NULL or character. For box containers. The value of \code{fill} (not
##' always respected) is used to control if expansion happens
##' vertically (\code{y}), horizontally (\code{x}) or both
##' (\code{both} or \code{TRUE}). For vertically filled box
##' containers, children always fill horizontally (atleast) and for
##' horizontally filled box containers, children always fill
##' vertically (atleast). This is important to realize when trying to
##' size buttons, say.
##' @param anchor NULL or integer. For box containers. The anchor argument is used to
##' position the child within the parent when there is more space
##' allocated than the child requests. This is specified with a
##' Cartesian pair in {-1,0,1} x {-1, 0, 1}. 
##' @param ... passed on to the 
##' @export
##' @rdname add
add <- function(obj, child, expand=FALSE, fill=NULL, anchor=NULL, ...) UseMethod("add")

##' Basic S3 method for add
##'
##' @export
##' @rdname add
##' @S3method add default
##' @method add default
add.default <- function(obj, child, expand=FALSE, fill=NULL, anchor=NULL, ...) {
  if(!isExtant(obj))  return()

  ## second dispatch based on type of child
  .add <- function(child, obj, ...) UseMethod(".add")  
  .add.GMenu <- function(child, obj, ...) {
    stop("Parent must be gwindow instance to add a menu")
  }
  .add.GToolBar <- function(child, obj, ...) {
    stop("Parent must be gwindow instance to add a toolbar")
  }
  .add.GStatusbar <- function(child, obj, ...) {
    stop("Parent must be gwindow instance to add a statusbar")
  }
  .add.default <- function(child, obj, expand, fill, anchor, ...) obj$add_child(child, expand=expand, fill=fill, anchor=anchor, ...)

  .add(child, obj, expand=expand, fill=fill, anchor=anchor, ...)
  
}



##' Delete child object from parent
##'
##' Delete may or may not remove a child. This is toolkit
##' specific. It may also be tied up with garbage collection. To avoid
##' that, keep a reference to the child object before deleting.
##' @export
##' @rdname add
delete <- function(obj, child) UseMethod("delete")

##' Basic S3 method for add
##'
##' @export
##' @rdname add
##' @S3method delete GContainer
##' @method delete GContainer
delete.GContainer <- function(obj, child) {
  if(isExtant(obj))
    obj$remove_child(child)
}

##' Dispose of object
##'
##' Dispose of object, primarily a window though this is modified in
##' \code{GNoteBook} and \code{GText}.
##' @param obj object to dispose
##' @param ... passed along
##' @export
##' @rdname dispose
dispose <- function(obj, ...) UseMethod("dispose")

##' main dispose method. Calls dispose for GWindow
##'
##' @export
##' @rdname dispose
##' @method dispose GComponent
##' @S3method dispose GComponent
dispose.GComponent <- function(obj, ...) {
  if(isExtant(obj))
    dispose(getTopLevel(obj))
}

## XXX dispose.GNotebook removes page

##' Get underlying toolkit widget
##'
##' At times a user may wish to access the underlying toolkit
##' widget. Although this is not cross-platform, one often has access
##' to many more methods of the object, than through those provided by
##' gWidgets.
##' @param obj object
##' @export
##' @rdname getToolkitWidget
getToolkitWidget <- function(obj) UseMethod("getToolkitWidget")

##' Basic S3 method 
##'
##' @export
##' @rdname getToolkitWidget
##' @method getToolkitWidget default
##' @S3method getToolkitWidget default
getToolkitWidget.default <- function(obj) getWidget(obj)

##' Get underlying toolkit widget from widget slot. Used internally
##'
##' @export
##' @rdname getToolkitWidget
getWidget <- function(obj) UseMethod("getWidget")



##' method for getWidget
##'
##' @rdname getToolkitWidget
##' @export
##' @method getWidget GComponent
##' @S3method getWidget GComponent
getWidget.GComponent <- function(obj) getWidget(obj$widget)

## implement getWidget.RGtkObject <- function(obj) obj say

##' Get underlying toolkit object from block slot
##'
##' @rdname getToolkitWidget
##' @export
getBlock <- function(obj) UseMethod("getBlock")

##' S3 method for getBlock generic
##'
##' @rdname getToolkitWidget
##' @export
##' @method getBlock GComponent
##' @S3method getBlock GComponent
getBlock.GComponent <- function(obj) getBlock(obj$block)


##' S3 method for getBlock generic
##'
##' For GWindow, the block is NULL
##' @rdname getToolkitWidget
##' @export
##' @method getBlock GWindow
##' @S3method getBlock GWindow
getBlock.GWindow <- function(obj) obj$widget
             
##' Get toplevel window containing object
##'
##' @export
##' @rdname getToolkitWidget
getTopLevel <- function(obj) UseMethod("getTopLevel")

##' getTopLevel method for components
##'
##' @export
##' @rdname getToolkitWidget
##' @method getTopLevel GComponent
##' @S3method getTopLevel GComponent
getTopLevel.GComponent <- function(obj) {
  if(!is(obj, "GComponent"))
    stop("Must call getTopLevel with a GComponent object")
  if(is(obj$parent, "uninitializedField") || is.null(obj$parent))
    return(obj)
  else
    getTopLevel(obj$parent)
}


## Container methods

##' Add a spring to box containers
##'
##' A spring will separate the children packed in the box container
##' prior to the spring be added and those being added, pushing the
##' two as far apart as the allocated space will allow.
##' @return NULL
##' @export
##' @rdname methods
addSpring <- function(obj) UseMethod("addSpring")

##' basic S3 generic to dispatch on
##'
##' Add spring to GContainer class
##' @export
##' @rdname methods
##' @method addSpring GContainer
##' @S3method addSpring GContainer
addSpring.GContainer <- function(obj) {
  obj$add_spring()
}


##' Add a space to a box container objects
##'
##' Inserts a specific amount of space between the previously packed
##' child and next one.
##' @param obj GContainer object
##' @param value space in pixels to add
##' @return NULL
##' @export
##' @rdname methods
addSpace <- function(obj, value) UseMethod("addSpace")

##' basic S3 generic to dispatch on
##'
##' Add space to GContainer class
##' @export
##' @rdname methods
##' @method addSpace GContainer
##' @S3method addSpace GContainer
addSpace.GContainer <- function(obj, value) {
  obj$add_space(value)
}

