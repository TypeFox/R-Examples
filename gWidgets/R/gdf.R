##' @include gtable.R

##' Class for a data frame editor
setClass("gDf",
         contains="gGridComponent",
         prototype=prototype(new("gGridComponent"))
         )

##' Constructor for a data frame editor
##'
##' @exports
##' @param items data frame to edit
##' @param name name of data frame (for saving, but not really needed)
##' @param do.subset logical. If \code{TRUE} a combobox to subset values by is produced.
##' @param container parent container
##' @param ... passed to container's \code{add} method
##' @param toolkit toolkit
##' @return An object of class \code{gDf} with overridden methods:
##'
##' \enumerate{
##'
##' \item \code{svalue} to return the selected cells
##'
##' \item \code{svalue<-} to set the selected cells
##'
##' \item \code{[} to get the current values
##'
##' \item \code{[<-} to set the current values. For most toolkits one
##' can not coerce the underlying class of a column (say from numeric
##' to character)
##'
##' \item \code{dim} dimension of object
##'
##' \item \code{length} number of columns
##'
##' }
##'
##' The widget may have underlying handlers assigned by \code{addHandlerClicked}
##'
##' }
##' 
gdf <- function(
                items = NULL, name = deparse(substitute(items)), do.subset = FALSE,
                container = NULL, ... ,
                toolkit=guiToolkit()){
  widget <- .gdf (toolkit,
                  items=items, name=name, do.subset=do.subset, container=container ,...
                  )
  obj <- new( 'gDf',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gdf
setGeneric( '.gdf' ,
           function(toolkit,
                    items = NULL, name = deparse(substitute(items)),
                    do.subset = FALSE,      container = NULL, ... )
           standardGeneric( '.gdf' ))


## methods
## addHandlerChanged, addHandlerClicked, addHandlerDoubleClicked
## addHandlerColumnClicked, addHandlerColumnDoubleClicked
