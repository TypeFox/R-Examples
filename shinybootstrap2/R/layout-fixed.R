#' Create a page with a fixed layout
#'
#' Functions for creating fixed page layouts. A fixed page layout consists of
#' rows which in turn include columns. Rows exist for the purpose of making sure
#' their elements appear on the same line (if the browser has adequate width).
#' Columns exist for the purpose of defining how much horizontal space within a
#' 12-unit wide grid it's elements should occupy. Fixed pages limit their width
#' to 940 pixels on a typical display, and 724px or 1170px on smaller and larger
#' displays respectively.
#'
#' @param ... Elements to include within the container
#' @param title The browser window title (defaults to the host URL of the page)
#' @param responsive \code{TRUE} to use responsive layout (automatically adapt
#'   and resize page elements based on the size of the viewing device)
#' @param theme Alternative Bootstrap stylesheet (normally a css file within the
#'   www directory). For example, to use the theme located at
#'   \code{www/bootstrap.css} you would use \code{theme = "bootstrap.css"}.
#'
#' @return A UI defintion that can be passed to the \link[shiny]{shinyUI} function.
#'
#' @details To create a fixed page use the \code{fixedPage} function and include
#'   instances of \code{fixedRow} and \code{\link{column}} within it. Note that
#'   unlike \code{\link{fluidPage}}, fixed pages cannot make use of higher-level
#'   layout functions like \code{sidebarLayout}, rather, all layout must be done
#'   with \code{fixedRow} and \code{column}.
#'
#' @note See the
#'   \href{https://github.com/rstudio/shiny/wiki/Shiny-Application-Layout-Guide}{
#'   Shiny Application Layout Guide} for additional details on laying out fixed
#'   pages.
#'
#' @seealso \code{\link{column}}
#'
#' @examples
#' library(shiny)
#'
#' fixedPage(
#'   title = "Hello, Shiny!",
#'   fixedRow(
#'     column(width = 4,
#'       "4"
#'     ),
#'     column(width = 3, offset = 2,
#'       "3 offset 2"
#'     )
#'   )
#' )
#'
#' @rdname fixedPage
#' @export
fixedPage <- function(..., title = NULL, responsive = TRUE, theme = NULL) {
  bootstrapPage(div(class = "container", ...),
                title = title,
                responsive = responsive,
                theme = theme)
}

#' @rdname fixedPage
#' @export
fixedRow <- function(...) {
  div(class = "row", ...)
}
