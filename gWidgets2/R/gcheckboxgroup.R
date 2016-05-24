##' @include methods.R
NULL

##' Constructor for checkbox group. A linked group of checkboxes, but not exclusive.
##'
##' @export
##' @param items checkbox labels
##' @param checked logical. Are values checked
##' @param horizontal logical. If true displayed horizontally, else vertically
##' @param use.table logical. If supported, and \code{TRUE} then uses a table widget with scrollbars
##' @inheritParams gwidget
##' @return Returns an object of class \code{GCheckboxGroup} for which
##' the following methods are overridden:
##' 
##' \itemize{
##' 
##' \item{ \code{svalue} Return the selected values or an empty
##' character vector. If \code{index=TRUE}, returns indices of
##' selected values.}
##' 
##' \item{ \code{svalue<-} Set the selected values one of three ways:
##' by label name, by a logical variable indicating which are selected
##' (if ambigous, logical wins), if \code{index=TRUE} by the indices
##' to select.}
##' 
##' \item{ \code{[} returns labels}
##' 
##' \item{ \code{[<-} set the label values. Should be able to shorten
##' or lengthen list}
##' 
##' }
##' @example inst/examples/ex-selectionwidgets.R
gcheckboxgroup <- function(
                           items, checked = FALSE, horizontal = FALSE,
                           use.table=FALSE, handler = NULL,
                           action = NULL, container = NULL, ... ,
                           toolkit=guiToolkit()){

  if(missing(items))
    items <- character(0)
  horizontal <- as.logical(horizontal)
  
  
  obj <- .gcheckboxgroup (toolkit,
                          items=items, checked=checked, horizontal=horizontal, use.table=use.table,
                          handler=handler, action=action, container=container, ...
                          )

  check_return_class(obj,  c("GCheckboxGroup","GCheckboxGroupTable"))
  obj
  
 
}


##' generic for toolkit dispatch
##' 
##' @rdname gcheckboxgroup
##' @export
.gcheckboxgroup <- function(toolkit,
                            items, checked = FALSE, horizontal = FALSE, use.table=FALSE,
                            handler = NULL, action = NULL,
                            container = NULL, ... ) UseMethod( '.gcheckboxgroup' )





##' Change handler for \code{GCheckboxGroup}
##'
##' Change handler for a \code{GCheckboxGroup}p is called when any of the
##' checkboxes changes state.
##'
##' @param obj receiver object
##' @export
##' @rdname gcheckboxgroup
##' @method addHandlerChanged GCheckboxGroup
##' @S3method addHandlerChanged GCheckboxGroup
addHandlerChanged.GCheckboxGroup <- function(obj, handler, action=NULL, ...) NextMethod()



##' svalue method
##'
##' The \code{svalue} methods refer to the selected values. By default
##' these are the item values, coerced to characterq. When
##' \code{index=TRUE} is specified, then the index is returned as an
##' integer vector. For setting, one may also use a vector of logicals
##' (which is recycled) for the index.
##'
##' @inheritParams svalue
##' @export
##' @rdname gcheckboxgroup
##' @method svalue GCheckboxGroup
##' @S3method svalue GCheckboxGroup
"svalue.GCheckboxGroup" <- function(obj, index=NULL,  drop=NULL, ...)   NextMethod()
