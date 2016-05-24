##' @include ws-model.R
NULL


##' Constructor for workspace variable browser
##'
##' A workspace browser widget.  The workspace browser displays values
##' in the global environment.  Displayed objects are shown in
##' categories.
##'
##' For defining the categories, the reference method \code{set_filter_classes}
##' takes a named list, the names defining the categories, the values
##' being the classes belonging to that category. Non categorized
##' values appear separately.  The default is defined in
##' \code{gWidgets2:::gvarbrowser_default_classes}.
##'
##' The variable browser uses an instance of \code{WSWatcherModel} to
##' monitor the global workspace. This instance may be useful for
##' other purposes. (For example, one may add an observer that is
##' called to listen for changes to the set of available data
##' frames.). The instance is available through the \code{ws_model}
##' property.
##'
##' The \code{svalue} method returns the selected variable names. If
##' \code{drop=FALSE} is given, the objects are returned.
##'
##' The widget should support dragging from without needing to specify
##' a \code{drag_source}, though this may be overridden.
##' 
##' Use \code{addHandlerChanged} to listen to activation of a variable
##' (double clicking). Use \code{addHandlerSelectionChanged} to
##' monitor change of selection.
##' @inheritParams gwidget
##' @export
##' @rdname gvarbrowser
gvarbrowser <- function(
                        handler = NULL,
                        action = "summary",
                        container = NULL ,...,
                        toolkit=guiToolkit()){
  obj <- .gvarbrowser (toolkit,
                          handler=handler, action=action, container = container, ...
                          )
  check_return_class(obj, "GVarBrowser")
  return(obj)
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gvarbrowser
.gvarbrowser <-  function(toolkit,
                          handler = NULL,
                          action = "summary", container = NULL,... )
  UseMethod( '.gvarbrowser' )



## Can override via option gWidgets2:gvarbrowser_classes, or set filter_classes
gvarbrowser_default_classes <-
  list("Data"=c("integer", "numeric",  "matrix", "character", "factor", "logical"),
       "Data sets"=c("data.frame", "list"),
       "Models"=c("lm", "rlm"),
       "Functions"=c("function")
       )


##' svalue method
##'
##' Return selected objects a string (when \code{drop=TRUE}) with
##' recursive values separated by \code{$}, or the objects themselves (when \code{drop=FALSE}).
##' @inheritParams svalue
##' @export
##' @rdname gvarbrowser
##' @method svalue GVarBrowser
##' @S3method svalue GVarBrowser
svalue.GVarBrowser <-  function(obj, index=FALSE, drop=TRUE, ...) NextMethod()


##################################################

##' Get an object from an environment specified by a string.
##'
##' @param value A single character value dispatches to get. For a
##' length 2 or more, then assumes object is recursive and extracts
##' named components
##' @param envir environment to look for values.
##' @return the object or an error
get_object_from_string <- function(value, envir=.GlobalEnv) {
  val <- get(value[1], envir)
  if(length(value) > 1) {
    f <- function(x, key) x[[key]]
    val <- Reduce(f, list(val, value[-1]))
  }
  val
}


##' Provide a short summary for an object
##'
##' @param x object
##' @export
##' @rdname short_summary
short_summary <- function(x) UseMethod("short_summary")

##' method for generic
##'
##' @inheritParams short_summary
##' @export
##' @rdname short_summary
##' @method short_summary default
##' @S3method short_summary default
short_summary.default <- function(x) sprintf("Object with class %s", class(x)[1])

##' method for generic
##'
##' @inheritParams short_summary
##' @export
##' @rdname short_summary
##' @method short_summary numeric
##' @S3method short_summary numeric
short_summary.numeric <- function(x) sprintf("Numeric object, length %s", length(x))

##' method for generic
##'
##' @inheritParams short_summary
##' @export
##' @rdname short_summary
##' @method short_summary character
##' @S3method short_summary character
short_summary.character <- function(x) sprintf("Character object, length %s", length(x))

##' method for generic
##'
##' @inheritParams short_summary
##' @export
##' @rdname short_summary
##' @method short_summary logical
##' @S3method short_summary logical
short_summary.logical <- function(x) sprintf("Logical object, length %s", length(x))

##' method for generic
##'
##' @inheritParams short_summary
##' @export
##' @rdname short_summary
##' @method short_summary data.frame
##' @S3method short_summary data.frame
short_summary.data.frame <- function(x) sprintf("Data frame, %s variables %s cases", ncol(x), nrow(x))

##' method for generic
##'
##' @inheritParams short_summary
##' @export
##' @rdname short_summary
##' @method short_summary matrix
##' @S3method short_summary matrix
short_summary.matrix <- function(x) sprintf("Matrix, %s by %s", ncol(x), nrow(x))

##' method for generic
##'
##' @inheritParams short_summary
##' @export
##' @rdname short_summary
##' @method short_summary list
##' @S3method short_summary list
short_summary.list <- function(x) sprintf("List with %s components", length(x))

##' method for generic
##'
##' @inheritParams short_summary
##' @export
##' @rdname short_summary
##' @method short_summary lm
##' @S3method short_summary lm
short_summary.lm <- function(x) sprintf("Model object")

##' method for generic
##'
##' @inheritParams short_summary
##' @export
##' @rdname short_summary
##' @method short_summary function
##' @S3method short_summary function
short_summary.function <- function(x) "Function"
