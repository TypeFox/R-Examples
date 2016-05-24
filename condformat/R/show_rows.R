#' Selects the rows to be printed
#'
#' Keeps the rows you mention in the printed table.
#' Compared to \code{\link[dplyr]{filter}}, show_rows does not remove the
#' rows from the actual data frame, they are removed only for printing.
#'
#' @param ... Comma separated list of unquoted extensions.
#'
#' @return A condformat_show_rows object, usually to be added to a condformat_tbl object
#'         as shown in the examples
#' @examples
#' library(condformat)
#' data(iris)
#' x <- head(iris)
#' condformat(x) + show_rows(Sepal.Length > 4.5, Species == "setosa")
#' @importFrom lazyeval lazy_dots
#' @export
#' @seealso \code{\link[dplyr]{filter}}
show_rows <- function(...) {
#   if (missing(row_names)) {
#     row_names = NA
#   }
  show_rows_(.dots = lazyeval::lazy_dots(...))  #,row_names = row_names)
}

#' @rdname show_rows
#' @inheritParams dplyr::filter
#' @importFrom lazyeval all_dots
#' @aliases show_rows
#' @export
#' @examples
#' library(condformat)
#' data(iris)
#' x <- head(iris)
#' condformat(x) + show_rows_(.docts = c("Sepal.Length > 4.5", "Species == 'setosa'"))
show_rows_ <- function(..., .dots) {
  dots <- lazyeval::all_dots(.dots, ...)
#   if (missing(row_names)) {
#     row_names <- NA
#   }
  showobj <- structure(list(row_expr = dots),
                       class = c("condformat_show_rows", "condformat_show_rows_filter"))
  return(showobj)
}

#' @importFrom dplyr filter_
render_show.condformat_show_rows_filter <- function(showobj, finalshow, x, ...) {

  # col_to_show: The columns that this show_columns would keep:
  xfiltered <- dplyr::filter_(x, .dots = showobj$row_expr)

  finalshow$xfiltered <- xfiltered
  return(finalshow)
}
