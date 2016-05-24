#' Create a header panel
#'
#' Create a header panel containing an application title.
#'
#' @param title An application title to display
#' @param windowTitle The title that should be displayed by the browser window.
#'   Useful if \code{title} is not a string.
#' @return A headerPanel that can be passed to \link{pageWithSidebar}
#'
#' @examples
#' headerPanel("Hello Shiny!")
#' @export
headerPanel <- function(title, windowTitle=title) {
  tagList(
    tags$head(tags$title(windowTitle)),
    div(class="span12", style="padding: 10px 0px;",
      h1(title)
    )
  )
}

#' Create a panel containing an application title.
#'
#' @param title An application title to display
#' @param windowTitle The title that should be displayed by the browser window.
#'
#' @details Calling this function has the side effect of including a
#'   \code{title} tag within the head. You can also specify a page title
#'   explicitly using the `title` parameter of the top-level page function.
#'
#'
#' @examples
#' titlePanel("Hello Shiny!")
#'
#' @export
titlePanel <- function(title, windowTitle=title) {
  tagList(
    tags$head(tags$title(windowTitle)),
    h2(style = "padding: 10px 0px;", title)
  )
}

#' Create a sidebar panel
#'
#' Create a sidebar panel containing input controls that can in turn be passed
#' to \code{\link{sidebarLayout}}.
#'
#' @param ... UI elements to include on the sidebar
#' @param width The width of the sidebar. For fluid layouts this is out of 12
#'   total units; for fixed layouts it is out of whatever the width of the
#'   sidebar's parent column is.
#' @return A sidebar that can be passed to \code{\link{sidebarLayout}}
#'
#' @examples
#' # Sidebar with controls to select a dataset and specify
#' # the number of observations to view
#' sidebarPanel(
#'   selectInput("dataset", "Choose a dataset:",
#'               choices = c("rock", "pressure", "cars")),
#'
#'   numericInput("obs", "Observations:", 10)
#' )
#' @export
sidebarPanel <- function(..., width = 4) {
  div(class=paste0("span", width),
    tags$form(class="well",
      ...
    )
  )
}

#' Create a main panel
#'
#' Create a main panel containing output elements that can in turn be passed to
#' \code{\link{sidebarLayout}}.
#'
#' @param ... Output elements to include in the main panel
#' @param width The width of the main panel. For fluid layouts this is out of 12
#'   total units; for fixed layouts it is out of whatever the width of the main
#'   panel's parent column is.
#' @return A main panel that can be passed to \code{\link{sidebarLayout}}.
#'
#' @examples
#' library(shiny)
#' # Show the caption and plot of the requested variable against mpg
#' mainPanel(
#'    h3(textOutput("caption")),
#'    plotOutput("mpgPlot")
#' )
#' @export
mainPanel <- function(..., width = 8) {
  div(class=paste0("span", width),
    ...
  )
}
