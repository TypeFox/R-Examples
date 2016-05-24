##' @include guiComponents.R

##' Class for tabular data display: gdf, gtable 
setClass("gGridComponent",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )


##' svalue method for gdf or gtable widget
##'
##' Main property of a table is the selection
##' @param obj object
##' @param index if \code{index=NULL} or \code{FALSE} then values from
##' the table are given if \code{index=TRUE} then the indices are
##' returned. The value of \code{drop} deteremines what is returned.
##' @param drop if \code{NULL} or \code{TRUE} then when
##' \code{index=TRUE} the indices of the currently selected rows are
##' returned. If \code{index=FALSE} then the values from the chosen
##' column are returned. Otherwise, if \code{drop=FALSE} then when
##' \code{index=TRUE} then -- if the toolkit supports it -- a list with
##' components \code{rows} and \code{columns} describing rectangles
##' of selected values.
##' @return a vector or data frame of values, or a vector or list of
##' indices. If no selection returns \code{NULL}.
setMethod("svalue", signature(obj="gGridComponent"),
          function(obj, index=NULL, drop=NULL, ... ) {
            .svalue(obj@widget, obj@toolkit, ...,index=index, drop=drop)            
          })




##' set selection in gdf or gtable object
##'
##' @param obj
##' @param index if \code{NULL} or \code{FALSE} then a specification
##' by match on the chosen columns. If \code{TRUE} then a
##' specification of indices. This can be as a vector of row indices,
##' or as a list with row and column specifications.
##' @param value replacement value
setReplaceMethod("svalue", signature(obj="gGridComponent"),
          function(obj, index=NULL, ...,value) {
            .svalue(obj@widget, obj@toolkit, index=index, ...) <- value
            return(obj)
          })




## leftbracket -- column coercion an issue, resizing frames an issue, replacement method

## leftbracket<-

## dim, dimnames, length, names, ....

##################################################

##' class to display tabular data for selection
setClass("gTable",
         contains="gGridComponent",
         prototype=prototype(new("gGridComponent"))
         )

##' A constructor for displaying tabular data for selection
##'
##' @exports
gtable <- function(
                   items, multiple = FALSE, chosencol = 1, icon.FUN = NULL,
                   filter.column = NULL, filter.labels = NULL, filter.FUN = NULL,
                   handler = NULL, action = NULL, container = NULL, ... ,
                   toolkit=guiToolkit()){

  ## coerce items
  if(!missing(items)) {
    if (is.vector(items))
      items <- data.frame(Values=items, stringsAsFactors=FALSE)
    if(is.matrix(items))
      items <- data.frame(items, stringsAsFactors=FALSE)
  }
  
  widget <- .gtable (toolkit,
                     items=items, multiple=multiple, chosencol=chosencol,
                     icon.FUN=icon.FUN,
                     filter.column=filter.column, filter.labels=filter.labels, filter.FUN=filter.FUN,
                     handler=handler, action=action, container=container ,...
                     )
  obj <- new( 'gTable',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gtable
setGeneric( '.gtable' ,
           function(toolkit,
                    items, multiple = FALSE, chosencol = 1,
                    icon.FUN = NULL,
                    filter.column = NULL, filter.labels = NULL, filter.FUN = NULL,
                    handler = NULL, action = NULL, container = NULL, ... )
           standardGeneric( '.gtable' ))


## handler code in subclasses
## addHandlerChanged, addHandlerClicked
