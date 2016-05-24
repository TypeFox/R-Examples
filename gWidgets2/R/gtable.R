##' @include methods.R
NULL





##' A constructor for displaying tabular data for selection
##'
##' The tabular widget allows a user to select one (or more) row(s)
##' using the mouse or keyboard selection. The selected rows are the
##' main property and are returned by svalue through their key (from
##' the column specified by \code{chosen.col}), or by index. The
##' change handler changes on double-click event. Use add handler
##' click to respond to a change in selection.
##'
##' Many generic methods for data frames are implemented for
##' \code{gtable}. These include \code{[}, \code{[<-}, \code{length},
##' \code{names}, and \code{names<-}
##' @param items data.frame specifies items for selection. May be a vector, matrix or data frame
##' @param multiple logical allow  multiple selection
##' @param chosen.col which value from the row is returned by selection
##' @param icon.col NULL or integer. If latter, specifies column containing stock icon
##' @param tooltip.col NULL or integer. If latter, specifies column containing tooltip
##' @inheritParams gwidget
##' @return Returns an object of class \code{GTable}
##' @export
##' @examples
##' \dontrun{
##' w <- gwindow("gtable example", visible=FALSE)
##' g <- gvbox(cont=w)
##' tbl <- gtable(mtcars, cont=g, expand=TRUE, fill=TRUE)
##' addHandlerClicked(tbl, handler=function(h,...) sprintf("You selected %s", svalue(h$obj)))
##' visible(w) <- TRUE
##'
##' }
gtable <- function(
                   items,
                   multiple = FALSE,
                   chosen.col = 1,
                   icon.col = NULL,
                   tooltip.col=NULL,
                   handler = NULL, action = NULL,
                   container = NULL, ... ,
                   toolkit=guiToolkit()){


  deprecated_args=list(
    "chosencol"="Use chosen.col",
    "icon.FUN" = "Use icon.col to specify column of items storing stock icon name"
    )
  check_deprecated(deprecated_args)
  
  ## coerce items to data frame
  if(!missing(items)) {
    if (is.vector(items))
      items <- data.frame(Values=items, stringsAsFactors=FALSE)
    if(is.matrix(items))
      items <- data.frame(items, stringsAsFactors=FALSE)
  }

  obj <- .gtable (toolkit,
                  items=items,
                  multiple=multiple,
                  chosen.col=chosen.col,
                  icon.col = icon.col,
                  tooltip.col = tooltip.col,
                  handler=handler,
                  action=action,
                  container=container ,...
                     )
  check_return_class(obj, "GTable")
  return(obj)
}


##' .gtable generic for toolkit dispatch
##'
##' @export
##' @rdname gtable
.gtable <- function(toolkit,
                    items,
                    multiple = FALSE,
                    chosen.col = 1,
                    icon.col = NULL,
                    tooltip.col=NULL,
                    handler = NULL, action = NULL,
                    container = NULL,
                    ... )
           UseMethod( '.gtable' )





##' "svalue" method
##'
##' For gtable one can pass in row(s) to select by index (when index=TRUE) or by match among the values in the chosen column.
##' For setting by index, a value of 0L or integer(0) will clear the current selection
##' @inheritParams svalue
##' @export
##' @rdname gtable
##' @method svalue GTable
##' @S3method svalue GTable
svalue.GTable <- function(obj, index=NULL, ..., value) NextMethod()

##' "[" method
##'
##' For \code{GTable} objects the \code{[} and \code{[<-} methods are
##' (mostly) implemented. The \code{[} method allows one to access the
##' object using the regular matrix notation (but there is no list
##' notation, e.g. \code{$} or \code{[[}, defined). The \code{[<-}
##' method is available, but for most toolkits is restricted: one can
##' not add columns, add rows, remove columns, remove rows, or change
##' the type of the column. For all of these actions, one can reassign
##' the items being displayed through the idiom \code{obj[] <-
##' new_items}. This allows the widget to resize or redo the column
##' renderers.
##' @param x \code{GTable} object
##' @param i row index
##' @param j column index
##' @param drop do we drop when subsetting
##' @export
##' @rdname gtable
##' @method [ GTable
##' @S3method [ GTable
"[.GTable" <- function(x, i, j, ..., drop=TRUE) NextMethod()


##' addHandlerChanged method
##'
##' The change handler for \code{GTable} is called when the selection
##' changes. This is often the result of a click event (but need not
##' be), although we alias to \code{addHandlerClicked}. For double
##' click events, see \code{addHandlerDoubleclick}.
##' @inheritParams addHandler
##' @export
##' @rdname gtable
##' @method addHandlerChanged GTable
##' @S3method addHandlerChanged GTable
addHandlerChanged.GTable <- function(obj, handler, action=NULL, ...) NextMethod()


##' addHandlerDoubleclick method
##'
##' Double clicking is used to activate an item (single click is
##' selection). We also bind pressing the Return key on an item to
##' initiate this signal
##' @inheritParams addHandler
##' @export
##' @rdname gtable
##' @method addHandlerDoubleclick GTable
##' @S3method addHandlerDoubleclick GTable
addHandlerDoubleclick.GTable <- function(obj, handler, action=NULL, ...) NextMethod()

##' visible
##'
##' For \code{GTable}, visibility refers to which rows are currently
##' shown, not whether the widget itself is shown or hidden. (For the
##' latter, place the widget into a container and adjust that). One
##' can use this method to perform filtering by adjusting which rows
##' are displayed according to some criteria that returns a logical.
##' @inheritParams svalue
##' @export
##' @rdname gtable
##' @method visible GTable
##' @S3method visible GTable
visible.GTable <- function(obj, ...) NextMethod()

##' Set size of table widget
##'
##' For \code{GTable} the \code{size<-}  method is overridden to allow
##' one to specify the column widths. To do so, pass in the values for
##' \code{width}, \code{height} or \code{column.widths} as named components
##' of a list. The value of \code{column.widths} should be a numeric
##' vector of pixel widths of length matching the number of columns.
##' @inheritParams size
##' @export
##' @usage \method{size}{GTable} (obj) <- value
##' @rdname gtable
##' @method size<- GTable
##' @S3method size<- GTable
"size<-.GTable" <- function(obj, value) NextMethod()
