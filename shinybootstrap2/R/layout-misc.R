#' Create a column within a UI definition
#'
#' Create a column for use within a  \code{\link{fluidRow}} or
#' \code{\link[shiny]{fixedRow}}
#'
#' @param width The grid width of the column (must be between 1 and 12)
#' @param ... Elements to include within the column
#' @param offset The number of columns to offset this column from the end of the
#'   previous column.
#'
#' @return A column that can be included within a
#'   \code{\link{fluidRow}} or \code{\link[shiny]{fixedRow}}.
#'
#'
#' @seealso \code{\link{fluidRow}}, \code{\link[shiny]{fixedRow}}.
#'
#' @examples
#' library(shiny)
#'
#' fluidRow(
#'   column(4,
#'     sliderInput("obs", "Number of observations:",
#'                 min = 1, max = 1000, value = 500)
#'   ),
#'   column(8,
#'     plotOutput("distPlot")
#'   )
#' )
#'
#' fluidRow(
#'   column(width = 4,
#'     "4"
#'   ),
#'   column(width = 3, offset = 2,
#'     "3 offset 2"
#'   )
#' )
#' @export
column <- function(width, ..., offset = 0) {

  if (!is.numeric(width) || (width < 1) || (width > 12))
    stop("column width must be between 1 and 12")

  colClass <- paste0("span", width)
  if (offset > 0)
    colClass <- paste0(colClass, " offset", offset)
  div(class = colClass, ...)
}


#' Layout a sidebar and main area
#'
#' Create a layout with a sidebar and main area. The sidebar is displayed with a
#' distinct background color and typically contains input controls. The main
#' area occupies 2/3 of the horizontal width and typically contains outputs.
#'
#' @param sidebarPanel The \link{sidebarPanel} containing input controls
#' @param mainPanel The \link{mainPanel} containing outputs
#' @param position The position of the sidebar relative to the main area ("left"
#'   or "right")
#' @param fluid \code{TRUE} to use fluid layout; \code{FALSE} to use fixed
#'   layout.
#'
#' @examples
#' library(shiny)
#'
#' # Define UI
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
#' @export
sidebarLayout <- function(sidebarPanel,
                          mainPanel,
                          position = c("left", "right"),
                          fluid = TRUE) {

  # determine the order
  position <- match.arg(position)
  if (position == "left") {
    firstPanel <- sidebarPanel
    secondPanel <- mainPanel
  }
  else if (position == "right") {
    firstPanel <- mainPanel
    secondPanel <- sidebarPanel
  }

  # return as as row
  if (fluid)
    fluidRow(firstPanel, secondPanel)
  else
    fixedRow(firstPanel, secondPanel)
}

#' Lay out UI elements vertically
#'
#' Create a container that includes one or more rows of content (each element
#' passed to the container will appear on it's own line in the UI)
#'
#' @param ... Elements to include within the container
#' @param fluid \code{TRUE} to use fluid layout; \code{FALSE} to use fixed
#'   layout.
#'
#' @seealso \code{\link{fluidPage}}, \code{\link[shiny]{flowLayout}}
#'
#' @examples
#' library(shiny)
#'
#' fluidPage(
#'   verticalLayout(
#'     a(href="http://example.com/link1", "Link One"),
#'     a(href="http://example.com/link2", "Link Two"),
#'     a(href="http://example.com/link3", "Link Three")
#'   )
#' )
#' @export
verticalLayout <- function(..., fluid = TRUE) {
  lapply(list(...), function(row) {
    col <- column(12, row)
    if (fluid)
      fluidRow(col)
    else
      fixedRow(col)
  })
}
