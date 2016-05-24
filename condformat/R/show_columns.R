#' Selects the variables to be printed
#'
#' Keeps the variables you mention in the printed table.
#' Compared to \code{\link[dplyr]{select}}, show_columns does not remove the
#' columns from the data frame, so formatting rules can still depend
#' on them.
#' @param ... Comma separated list of unquoted extensions.
#' @param col_names Character vector with the column names for the selected columns
#'
#' @return A condformat_show_columns object, usually to be added to a condformat_tbl object
#' @examples
#' data(iris)
#' x <- head(iris)
#' condformat(x) + show_columns(Sepal.Length, Sepal.Width, Species)
#' condformat(x) + show_columns(-Petal.Length, -Petal.Width)
#' condformat(x) + show_columns(starts_with("Petal"), Species)
#' @importFrom lazyeval lazy_dots
#' @export
#' @seealso \code{\link[dplyr]{select}}
show_columns <- function(..., col_names) {
  if (missing(col_names)) {
    col_names = NA
  }
  show_columns_(.dots = lazyeval::lazy_dots(...), col_names = col_names)
}

#' @rdname show_columns
#' @inheritParams dplyr::select
#' @importFrom lazyeval all_dots
#' @export
#' @examples
#' data(iris)
#' x <- head(iris)
#' condformat(x) + show_columns_(.dots = c("Sepal.Length", "Species"))
#' condformat(x) + show_columns_(.dots = c("Sepal.Length", "Species"),
#'                               col_names = c("Sepal Length", "Species"))
show_columns_ <- function(..., .dots, col_names) {
  dots <- lazyeval::all_dots(.dots, ...)
  if (missing(col_names)) {
    col_names <- NA
  }
  showobj <- structure(list(column_expr = dots, col_names = col_names),
                       class = c("condformat_show_columns",
                                 "condformat_show_columns_select"))
  return(showobj)
}

#' @importFrom dplyr select_vars_
render_show.condformat_show_columns_select <- function(showobj, finalshow, x, ...) {

  # col_to_show: The columns that this show_columns would keep:
  col_to_show <- dplyr::select_vars_(colnames(x), showobj$column_expr)

  # Assign the names we want to use for those columns:
  if (!identical(showobj$col_names, NA)) {
    names(col_to_show) <- showobj$col_names
  } else {
    names(col_to_show) <- col_to_show
  }

  # If a variable had already been excluded, do not show it:
  col_to_show <- col_to_show[col_to_show %in% finalshow$cols]
  finalshow$cols <- col_to_show
  return(finalshow)
}
