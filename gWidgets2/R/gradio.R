##' @include methods.R
NULL


##' Constructor for radio button widget
##'
##' A radio button group allows a user to select one from many
##' items. In \pkg{gWidgets2} the radio button widget shows 2 or more
##' items. The items are coerced to characters, usually by the
##' underlying toolkit. Use the \code{coerce_with} property to set a
##' function, such as \code{as.numeric}, to coerce the return value
##' during the \code{svalue} code. The items are referred to with the
##' \code{[} method, the selected one with \code{svalue}.
##' @param items items to select from
##' @param selected index of initially selected item
##' @param horizontal layout direction
##' @inheritParams gwidget
##' @export
##' @rdname gradio
##' @example inst/examples/ex-selectionwidgets.R
gradio <- function(items,selected=1, horizontal=FALSE,
                   handler=NULL, action=NULL,
                   container=NULL, ...,
                   toolkit=guiToolkit()) {
  
  ## check input
  if(length(x <- unique(items) ) != length(items))
    message("Using unique items for selection values")
  
  obj <- .gradio(toolkit, x, selected, horizontal, handler, action, container,...)
  check_return_class(obj, "GRadio")
  obj
}

##' Generic for method dispatch
##'
##' @export
##' @rdname gradio
.gradio <- function(toolkit,
                    items, selected=1, horizontal=FALSE, handler=NULL, action=NULL,
                    container=NULL,
                    ...) UseMethod(".gradio")



##' svalue method
##'
##' The svalue method returns the radio button label or its index if
##' \code{index=TRUE}. Labels are coerced to character by many of the
##' toolkits. To be sure to return a numeric value, one can assign to
##' the \code{coerce_with} property, e.g., \code{obj$coerce_with <-
##' as.numeric}. For all widgets, if a function is specified to
##' \code{coerce_with}  it will be called on the value returned by
##' \code{svalue}.
##' @inheritParams svalue
##' @export
##' @rdname gradio
##' @S3method svalue GRadio
##' @method svalue GRadio
svalue.GRadio <- function(obj, index=NULL, drop=TRUE, ...) NextMethod()


##' svalue<- method
##'
##' For a radio button group, for \code{svalue} the value can be
##' referred to by index or label. 
##' @inheritParams svalue
##' @export
##' @usage \method{svalue}{GRadio} (obj,index=NULL,drop=TRUE,...) <- value
##' @rdname gradio
##' @S3method svalue<- GRadio
##' @method svalue<- GRadio
"svalue<-.GRadio" <- function(obj, index=NULL, drop=TRUE, ..., value) {
  if(!is.null(index) && index) {
    value <- as.integer(value)[1]
    if(value < 1 || value > length(obj)) warning(gettext("Index is out of range"))
  }
  if(is.null(index) || !index) {
    if(! value %in% obj[])
      warning(gettext("Value specified is not one of the items"))
  }

  NextMethod()
}


##' assign items for gradio
##'
##' Check for repeated items before passing on to \code{set_items}
##' @param x \code{GRadio} object
##' @param i button index. Leavel as missing to replace items to select from.
##' @param j ignored
##' @param value items to assigns a choices for the buttons
##' @export
##' @usage \method{[}{GRadio} (x, i, j, ...) <- value
##' @rdname gradio
##' @method [<- GRadio
##' @S3method [<- GRadio
"[<-.GRadio" <- function(x, i, j, ..., value) {
  ## check input
  if(length(value) != length(value <- unique(value)))
    message("Using unique values for selection values")
  NextMethod()
}
