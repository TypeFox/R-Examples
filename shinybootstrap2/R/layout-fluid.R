#' Create a page with fluid layout
#'
#' Functions for creating fluid page layouts. A fluid page layout consists of
#' rows which in turn include columns. Rows exist for the purpose of making sure
#' their elements appear on the same line (if the browser has adequate width).
#' Columns exist for the purpose of defining how much horizontal space within a
#' 12-unit wide grid it's elements should occupy. Fluid pages scale their
#' components in realtime to fill all available browser width.
#'
#' @param ... Elements to include within the page
#' @param title The browser window title (defaults to the host URL of the page).
#'   Can also be set as a side effect of the \code{\link{titlePanel}} function.
#' @param responsive \code{TRUE} to use responsive layout (automatically adapt
#'   and resize page elements based on the size of the viewing device)
#' @param theme Alternative Bootstrap stylesheet (normally a css file within the
#'   www directory). For example, to use the theme located at
#'   \code{www/bootstrap.css} you would use \code{theme = "bootstrap.css"}.
#'
#' @return A UI defintion that can be passed to the \link[shiny]{shinyUI} function.
#'
#' @details To create a fluid page use the \code{fluidPage} function and include
#'   instances of \code{fluidRow} and \code{\link{column}} within it. As an
#'   alternative to low-level row and column functions you can also use
#'   higher-level layout functions like \code{\link{sidebarLayout}}.
#'
#' @note See the
#'   \href{https://github.com/rstudio/shiny/wiki/Shiny-Application-Layout-Guide}{
#'    Shiny-Application-Layout-Guide} for additional details on laying out fluid
#'   pages.
#'
#' @seealso \code{\link{column}}, \code{\link{sidebarLayout}}
#'
#' @examples
#' library(shiny)
#'
#' fluidPage(
#'
#'   # Application title
#'   titlePanel("Hello Shiny!"),
#'
#'   sidebarLayout(
#'
#'     # Sidebar with a slider input
#'     sidebarPanel(
#'       sliderInput("obs",
#'                   "Number of observations:",
#'                   min = 0,
#'                   max = 1000,
#'                   value = 500)
#'     ),
#'
#'     # Show a plot of the generated distribution
#'     mainPanel(
#'       plotOutput("distPlot")
#'     )
#'   )
#' )
#'
#' fluidPage(
#'   title = "Hello Shiny!",
#'   fluidRow(
#'     column(width = 4,
#'       "4"
#'     ),
#'     column(width = 3, offset = 2,
#'       "3 offset 2"
#'     )
#'   )
#' )
#'
#' @rdname fluidPage
#' @export
fluidPage <- function(..., title = NULL, responsive = TRUE, theme = NULL) {
  bootstrapPage(div(class = "container-fluid", ...),
                title = title,
                responsive = responsive,
                theme = theme)
}


#' @rdname fluidPage
#' @export
fluidRow <- function(...) {
  div(class = "row-fluid", ...)
}
