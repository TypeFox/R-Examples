##' @include methods.R
NULL

##' constructor for a combobox
##'
##' A combobox can be either a drop down list (editable=FALSE), or a drop-down list and edit area (a combobox).
##' @param items Items to select from. A vector or a data frame. If a
##' data frame, then first column is values. Second is optional, but
##' can specify a stock icon name, third is optional and can be used
##' to specify a tooltip. These may not be supported in all toolkits.
##' @param selected integer. Which item (by index) is selected. Use -1 for no selection
##' @param editable logical. Is user allowed to edit value
##' @param coerce.with A function of function name to be called before
##' selected value is returned by \code{svalue}
##' @param handler Called when combobox value is changed.
##' @param action passed to handler
##' @param container parent container
##' @param ... passed to parent container's \code{add} method
##' @param toolkit toolkit
##' @return Returns an object of class \code{GComboBox} for which the following methods are overriden:
##' \enumerate{
##' \item \code{svalue} Return selected value by name or (if \code{index=TRUE} by index). The latter only if \code{editable=FALSE}.
##' 
##' \item \code{svalue<-} Set the selected value by value or if \code{index=TRUE} by index.
##'
##' \item \code{[} return items to select from
##'
##' \item \code{[<-} Set items to select from.
##' }
##' @export
gcombobox <- function(
                      items, selected = 1, editable = FALSE, coerce.with=NULL,
                      handler = NULL, action = NULL, container = NULL, ... ,
                      toolkit=guiToolkit()) {
  
  items <- .make_gcombobox_items(items)
  
  
  obj <- .gcombobox(toolkit,
                    items=items, selected=selected, editable=editable, coerce.with=coerce.with,
                    handler=handler, action=action, container=container, ...
                    )

  check_return_class(obj,"GComboBox")
  obj
  
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gcombobox
.gcombobox <- function(toolkit,
                       items, selected = 1, editable = FALSE, coerce.with = NULL,
                       handler = NULL,  action = NULL, container = NULL, ... ) UseMethod( '.gcombobox' )

##' Alias for gcombobox. Deprecated
##'
##' @export
##' @rdname gcombobox
gdroplist <- function(...) {
  .Deprecated("Use gcombobox, not gdroplist")
  gcombobox(...)
}




##' Helper function
##'
##' Non exported helper function to coerce items into a data frame.
##' First column contains the values, second stock icons, third tooltips
##' @rdname gcombobox
.make_gcombobox_items <- function(value) {
  if(!is.data.frame(value) && !is.matrix(value))
    value <- data.frame(value, stringsAsFactors=FALSE)
  
  if(!is.data.frame(value))
    value <- as.data.frame(value)

  value
}



##' replacement method for combobox selection items
##'
##' Ensure that value is a data frame. One can pass a vector or a
##' one-column data frame to inidicate the possible values for
##' selection, a second column is used for an icons (if possible), a
##' third for a tooltip (if possible).
##'
##' @param x combobox object
##' @param i item index
##' @param j ignored
##' @param value new items for selection
##' @export
##' @usage \method{[}{GComboBox} (x, i , j, ...) <- value
##' @rdname gcombobox
##' @method [<- GComboBox
##' @S3method [<- GComboBox
"[<-.GComboBox" <- function(x, i, j, ..., value) {

  value <- .make_gcombobox_items(value)
  
  NextMethod()
}
  


##' Change handler for \code{GcomboBox}
##'
##' Change handler for a non-editabe combobox is called when a new
##' value is selected. For editable comboboxes, the handler is also
##' called when the text entry box is activated.
##'
##' @inheritParams addHandlerChanged
##' @export
##' @rdname gcombobox
##' @method addHandlerChanged GComboBox
##' @S3method addHandlerChanged GComboBox
addHandlerChanged.GComboBox <- function(obj, handler, action=NULL, ...) NextMethod()



##' svalue method
##'
##' The \code{svalue} method for a combobox object refers to its main
##' property, the selected value. When \code{index=FALSE} (or
##' \code{NULL}) the value is returned. If \code{index=TRUE} the index
##' of the object within the set of items is used.
##' @inheritParams svalue
##' @export
##' @rdname gcombobox
##' @method svalue GComboBox
##' @S3method svalue GComboBox
svalue.GComboBox <- function(obj, index=NULL, drop=NULL, ...)   NextMethod()

