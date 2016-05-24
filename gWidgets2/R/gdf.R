##' @include methods.R
NULL

##' Constructor for a data frame editor
##'
##' Implementation varies wildly, but should provide at minimum
##' functionality of \code{edit.data.frame}. A single mouse click on a
##' cell should select that cell, a double click should initiate
##' editing of that cell.
##'
##' Contents of the data frame can be accessed via \code{[} and manipulated with \code{[<-}.
##' 
##' The \code{save_data} reference class method can be called to
##' save the data into a variable with the specified name.
##'
##' example in \code{inst/examples/ex-gdf.R}
##' @param items data frame to edit
##' @param container parent container
##' @param handler called on cell change
##' @param action passed to handler
##' @param ... passed to container's \code{add} method
##' @param toolkit toolkit
##' @export
##' @return An object of class \code{gDf}.
gdf <- function(
    
                items = NULL, 
                handler=NULL, action=NULL,
                container = NULL, ... ,
                toolkit=guiToolkit()){
  obj <- .gdf (toolkit,
               items=items, 
               handler=handler, action=action, container=container ,...
               )

  deprecated_args <- list(do.subset=c("do.subset argument has been deprecated.","See examples for how to create that functionality"),
                          name=c("name argument is deprecated, as this should be managed by programmer")
                          )


  check_deprecated(deprecated_args, ...)

  check_return_class(obj, "GDf")
  return(obj)
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gdf
.gdf <- function(toolkit,
                 items = NULL,
                 handler=NULL, action=NULL, container = NULL, ... )
  UseMethod( '.gdf' )


## methods
## addHandlerChanged, addHandlerClicked, addHandlerDoubleClicked
## addHandlerColumnClicked, addHandlerColumnDoubleClicked

##' Change handler for GDf
##'
##' Assign handler to be called when a cell, row or column changes
##' @inheritParams addHandler
##' @export
##' @rdname gdf
##' @method addHandlerChanged GDf
##' @S3method addHandlerChanged GDf
addHandlerChanged.GDf <- function(obj, handler, action=NULL, ...) NextMethod()
  

##' "svalue" method
##'
##' For gdf svalue refers to the selected values.
##' @inheritParams svalue
##' @export
##' @rdname gdf
##' @method svalue GDf
##' @S3method svalue GDf
svalue.GDf <- function(obj, index=NULL, drop=TRUE, ...) NextMethod()


##' "visible" method
##'
##' \code{visible} is used to refer to which rows are being shown.
##' @inheritParams visible
##' @export
##' @usage \method{visible}{GDf} (obj) <- value
##' @rdname gdf
##' @method visible<- GDf
##' @S3method visible<- GDf
"visible<-.GDf" <- function(obj, value) NextMethod()


##' get current items
##'
##' The current items for a \code{gdf} object are both the visible and
##' non-visible items. To retrieve just the currently visible items,
##' use the idiom \code{obj[visible(obj), ]}.
##' @inheritParams gdf
##' @export
##' @rdname gWidgets2-S3methods
##' @method [ GDf
##' @S3method [ GDf
"[.GDf" <- function(x, i, j, ..., drop=TRUE) NextMethod()
