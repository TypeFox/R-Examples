#' Converts a data frame or matrix to a condformat_tbl class.
#'
#' condformat_tbl objects allow to add conditional formatting information,
#' that can be viewed when the condformat_tbl object is printed.
#'
#' They also have all properties from conventional data frames.
#'
#' @param x A matrix or data.frame
#'
#' @return The condformat_tbl object. This object can be used
#'         in a ggplot-like syntax to apply conditional formatting
#'         rules. It can also be used as a conventional
#'         data frame.
#'
#'         The condformat_tbl print method generates an htmlTable, to be
#'         viewed using RStudio Viewer or an HTML browser, as available.
#' @importFrom dplyr tbl_df
#' @examples
#' data(iris)
#' condformat(iris[1:5,])
#' @export
condformat <- function(x) {
  x <- dplyr::tbl_df(x)
  attr(x, "condformat") <- list(show = list(rows = list(),
                                            cols = list()),
                                rules = list(),
                                themes = list())

  class(x) <- c("condformat_tbl", class(x))
  return(x)
}


#' Combines the data with the formatting rules and graphical options.
#'
#' @param x A condformat_tbl object
#' @param obj A condformat_show or a condformat_rule object to be combined
#'            Any other type of object will be added as expected to the data frame.
#' @return x, with extended condformat_tbl attributes
#' @examples
#' data(iris)
#' condformat(iris) + show_columns(Species)
#' @method + condformat_tbl
#' @export
"+.condformat_tbl" <- function(x, obj) {
  if (inherits(obj, "condformat_show_columns")) {
    condformatopts <- attr(x, "condformat")
    condformatopts$show$cols <- c(condformatopts$show$cols, list(obj))
    attr(x, "condformat") <- condformatopts
    return(x)
  } else if (inherits(obj, "condformat_show_rows")) {
    condformatopts <- attr(x, "condformat")
    condformatopts$show$rows <- c(condformatopts$show$rows, list(obj))
    attr(x, "condformat") <- condformatopts
    return(x)
  } else if (inherits(obj, "condformat_rule")) {
    condformatopts <- attr(x, "condformat")
    condformatopts$rules <- c(condformatopts$rules, list(obj))
    attr(x, "condformat") <- condformatopts
    return(x)
  } else if (inherits(obj, "condformat_theme")) {
    condformatopts <- attr(x, "condformat")
    condformatopts$themes <- c(condformatopts$themes, list(obj))
    attr(x, "condformat") <- condformatopts
    return(x)
  } else {
    NextMethod()
  }
}
