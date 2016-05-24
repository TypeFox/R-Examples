##' @include guiComponentWithItems.R

##' checkboxgroup class
setClass("gCheckboxGroup",
         contains="guiComponentWithItems",
         prototype=prototype(new("guiComponentWithItems"))
         )

##' Constructor for checkbox group. A linked group of checkboxes, but not exclusive.
##'
##' @export
##' @param items checkbox labels
##' @param checked logical. Are values checked
##' @param horizontal logical. If true displayed horizontally, else vertically
##' @param use.table logical. If supported, and \code{TRUE} then use a table widget with scrollbars
##' @param handler Handler called when state toggles
##' @param action passed to handler when called
##' @param container parent container
##' @param ... passed to \code{add} method of parent
##' @param toolkit toolkit
##' @example ~/pmg/r-forge/gwidgets/pkg/gWidgets/inst/tests/
##' @return Returns an object of class \code{gCheckboxGroup} for which the following methods are overridden:
##' %
##' \example{
##'
##' \item \code{svalue} Return the selected values or an empty
##' character vector. If \code{index=TRUE}, returns indices of
##' selected values.
##' 
##' \item \code{svalue<-} Set the selected values one of three ways:
##' by label name, by a logical variable indicating which are selected
##' (if ambigous, logical wins), if \code{index=TRUE} by the indices
##' to select.
##' 
##' \item \code{[} returns labels
##' 
##' \item \code{[<-} set the label values. Should be able to shorten
##' or lengthen list
##' 
##' }
gcheckboxgroup <- function(
                           items, checked = FALSE, horizontal = FALSE,
                           use.table=FALSE, handler = NULL,
                           action = NULL, container = NULL, ... ,
                           toolkit=guiToolkit()){

  if(missing(items))
    items <- character(0)
  horizontal <- as.logical(horizontal)
  
  
  widget <- .gcheckboxgroup (toolkit,
    items=items, checked=checked, horizontal=horizontal, use.table=use.table,
    handler=handler, action=action, container=container, ...
    )
  obj <- new( 'gCheckboxGroup',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gcheckboxgroup
##' @export
setGeneric( '.gcheckboxgroup' ,
           function(toolkit,
                    items, checked = FALSE, horizontal = FALSE,
                    handler = NULL, action = NULL,
                    container = NULL, ... ) standardGeneric( '.gcheckboxgroup' )) 



