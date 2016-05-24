#' Create an icon
#'
#' Create an icon for use within a page. Icons can appear on their own,
#' inside of a button, or as an icon for a \code{\link{tabPanel}} within a
#' \code{\link{navbarPage}}.
#'
#' @param name Name of icon. Icons are drawn from the
#'   \href{http://fontawesome.io/icons/}{Font Awesome} library. Note that the
#'   "fa-" prefix should not be used in icon names (i.e. the "fa-calendar" icon
#'   should be referred to as "calendar")
#' @param class Additional classes to customize the style of the icon (see the
#'  \href{http://fontawesome.io/examples/}{usage examples} for
#'   details on supported styles).
#' @param lib Icon library to use ("font-awesome" is only one currently
#'   supported)
#'
#' @return An icon element
#'
#'
#' @examples
#' library(shiny)
#' icon("calendar")            # standard icon
#' icon("calendar", "fa-3x")   # 3x normal size
#'
#' # add an icon to a submit button
#' submitButton("Update View", icon = icon("refresh"))
#'
#' navbarPage("App Title",
#'   tabPanel("Plot", icon = icon("bar-chart-o")),
#'   tabPanel("Summary", icon = icon("list-alt")),
#'   tabPanel("Table", icon = icon("table"))
#' )
#'
#' @export
icon <- function(name, class = NULL, lib = "font-awesome") {

  # determine stylesheet
  if (!identical(lib, "font-awesome"))
    stop("Unknown font library '", lib, "' specified")

  # build the icon class (allow name to be null so that other functions
  # e.g. builtTabset can pass an explict class value)
  iconClass <- ""
  if (!is.null(name))
    iconClass <- paste0("fa fa-", name)
  if (!is.null(class))
    iconClass <- paste(iconClass, class)

  # return the element and css dependency
  tagList(
    singleton(tags$head(
      tags$link(rel = "stylesheet", type = "text/css",
                href = 'shared/font-awesome/css/font-awesome.min.css')
    )),
    tags$i(class = iconClass)
  )
}

# Helper funtion to extract the class from an icon
iconClass <- function(icon) {
  if (!is.null(icon)) icon[[2]]$attribs$class
}
