#' Create a Bootstrap page
#'
#' Create a Shiny UI page that loads the CSS and JavaScript for
#' \href{http://getbootstrap.com/2.3.2/}{Bootstrap}, and has no content in the
#' page body (other than what you provide).
#'
#' This function is primarily intended for users who are proficient in
#' HTML/CSS, and know how to lay out pages in Bootstrap. Most applications
#' should use \code{\link{fluidPage}} along with layout functions like
#' \code{\link{fluidRow}} and \code{\link{sidebarLayout}}.
#'
#' @param ... The contents of the document body.
#' @param title The browser window title (defaults to the host URL of the page)
#' @param responsive \code{TRUE} to use responsive layout (automatically adapt
#'   and resize page elements based on the size of the viewing device)
#' @param theme Alternative Bootstrap stylesheet (normally a css file within the
#'   www directory, e.g. \code{www/bootstrap.css})
#'
#' @return A UI defintion that can be passed to the \link[shiny]{shinyUI} function.
#'
#' @note The \code{basicPage} function is deprecated, you should use the
#'   \code{\link{fluidPage}} function instead.
#'
#' @seealso \code{\link{fluidPage}}, \code{\link{fixedPage}}
#'
#' @export
bootstrapPage <- function(..., title = NULL, responsive = TRUE, theme = NULL) {

  # required head tags for boostrap
  importBootstrap <- function(min = TRUE) {

    ext <- function(ext) {
      ifelse(min, paste(".min", ext, sep=""), ext)
    }
    cssExt <- ext(".css")
    jsExt = ext(".js")
    bs <- c(
      file = system.file("www/bootstrap", package = "shinybootstrap2")
    )

    list(
      htmlDependency("bootstrap", "2.3.2", bs,
        script = sprintf("js/bootstrap%s", jsExt),
        stylesheet = c(
          "shiny-bootstrap2.css",
          if (is.null(theme)) sprintf("css/bootstrap%s", cssExt)
        )
      ),
      if (responsive) {
        htmlDependency("bootstrap-responsive", "2.3.2", bs,
          stylesheet = sprintf("css/bootstrap-responsive%s", cssExt),
          meta = list(viewport = "width=device-width, initial-scale=1.0")
        )
      }
    )
  }

  attachDependencies(
    tagList(
      if (!is.null(title)) tags$head(tags$title(title)),
      if (!is.null(theme)) {
        tags$head(tags$link(rel="stylesheet", type="text/css", href = theme))
      },

      # remainder of tags passed to the function
      list(...)
    ),
    importBootstrap()
  )
}

#' @rdname bootstrapPage
#' @export
basicPage <- function(...) {
  bootstrapPage(div(class="container-fluid", list(...)))
}

#' Create a page with a sidebar
#'
#' Create a Shiny UI that contains a header with the application title, a
#' sidebar for input controls, and a main area for output.
#'
#' @param headerPanel The \link{headerPanel} with the application title
#' @param sidebarPanel The \link{sidebarPanel} containing input controls
#' @param mainPanel The \link{mainPanel} containing outputs

#' @return A UI defintion that can be passed to the \link[shiny]{shinyUI} function
#'
#' @note This function is deprecated. You should use \code{\link{fluidPage}}
#' along with \code{\link{sidebarLayout}} to implement a page with a sidebar.
#'
#' @examples
#' library(shiny)
#'
#' pageWithSidebar(
#'
#'   # Application title
#'   headerPanel("Hello Shiny!"),
#'
#'   # Sidebar with a slider input
#'   sidebarPanel(
#'     sliderInput("obs",
#'                 "Number of observations:",
#'                 min = 0,
#'                 max = 1000,
#'                 value = 500)
#'   ),
#'
#'   # Show a plot of the generated distribution
#'   mainPanel(
#'     plotOutput("distPlot")
#'   )
#' )
#'
#' @export
pageWithSidebar <- function(headerPanel,
                            sidebarPanel,
                            mainPanel) {

  bootstrapPage(
    # basic application container divs
    div(
      class="container-fluid",
      div(class="row-fluid",
          headerPanel
      ),
      div(class="row-fluid",
          sidebarPanel,
          mainPanel
      )
    )
  )
}

#' Create a page with a top level navigation bar
#'
#' Create a page that contains a top level navigation bar that can be used to
#' toggle a set of \code{\link[shiny]{tabPanel}} elements.
#'
#' @param title The title to display in the navbar
#' @param ... \code{\link[shiny]{tabPanel}} elements to include in the page
#' @param id If provided, you can use \code{input$}\emph{\code{id}} in your
#'   server logic to determine which of the current tabs is active. The value
#'   will correspond to the \code{value} argument that is passed to
#'   \code{\link[shiny]{tabPanel}}.
#' @param header Tag of list of tags to display as a common header above all
#'   tabPanels.
#' @param footer Tag or list of tags to display as a common footer below all
#'   tabPanels
#' @param inverse \code{TRUE} to use a dark background and light text for the
#'   navigation bar
#' @param collapsable \code{TRUE} to automatically collapse the navigation
#'   elements into a menu when the width of the browser is less than 940 pixels
#'   (useful for viewing on smaller touchscreen device)
#' @param fluid \code{TRUE} to use a fluid layout. \code{FALSE} to use a fixed
#'   layout.
#' @param responsive \code{TRUE} to use responsive layout (automatically adapt
#'   and resize page elements based on the size of the viewing device)
#' @param theme Alternative Bootstrap stylesheet (normally a css file within the
#'   www directory). For example, to use the theme located at
#'   \code{www/bootstrap.css} you would use \code{theme = "bootstrap.css"}.
#' @param windowTitle The title that should be displayed by the browser window.
#'   Useful if \code{title} is not a string.
#'
#' @return A UI defintion that can be passed to the \link[shiny]{shinyUI} function.
#'
#' @details The \code{navbarMenu} function can be used to create an embedded
#'   menu within the navbar that in turns includes additional tabPanels (see
#'   example below).
#'
#' @seealso \code{\link[shiny]{tabPanel}}, \code{\link[shiny]{tabsetPanel}}
#'
#' @examples
#' \dontrun{
#' navbarPage("App Title",
#'   tabPanel("Plot"),
#'   tabPanel("Summary"),
#'   tabPanel("Table")
#' )
#'
#' navbarPage("App Title",
#'   tabPanel("Plot"),
#'   navbarMenu("More",
#'     tabPanel("Summary"),
#'     tabPanel("Table")
#'   )
#' )
#' }
#' @export
navbarPage <- function(title,
                       ...,
                       id = NULL,
                       header = NULL,
                       footer = NULL,
                       inverse = FALSE,
                       collapsable = FALSE,
                       fluid = TRUE,
                       responsive = TRUE,
                       theme = NULL,
                       windowTitle = title) {

  # alias title so we can avoid conflicts w/ title in withTags
  pageTitle <- title

  # navbar class based on options
  navbarClass <- "navbar navbar-static-top"
  if (inverse)
    navbarClass <- paste(navbarClass, "navbar-inverse")

  # build the tabset
  tabs <- list(...)
  tabset <- buildTabset(tabs, "nav", NULL, id)

  # built the container div dynamically to support optional collapsability
  if (collapsable) {
    navId <- paste("navbar-", shiny:::p_randomInt(1000, 10000), sep="")
    containerDiv <- div(class="container",
                        tags$button(type="button",
                                    class="btn btn-navbar",
                                    `data-toggle`="collapse",
                                    `data-target`=".nav-collapse",
                          span(class="icon-bar"),
                          span(class="icon-bar"),
                          span(class="icon-bar")
                        ),
                        span(class="brand pull-left", pageTitle),
                        div(class="nav-collapse collapse",
                            id=navId,
                            tabset$navList),
                            tags$script(paste(
                              "$('#", navId, " a:not(.dropdown-toggle)').click(function (e) {
                                  e.preventDefault();
                                  $(this).tab('show');
                                  if ($('.navbar .btn-navbar').is(':visible'))
                                    $('.navbar .btn-navbar').click();
                               });", sep="")))
  } else {
    containerDiv <- div(class="container",
                        span(class="brand pull-left", pageTitle),
                        tabset$navList)
  }

  # create a default header if necessary
  if (length(header) == 0)
    header <- HTML("&nbsp;")

  # function to return plain or fluid class name
  className <- function(name) {
    if (fluid)
      paste(name, "-fluid", sep="")
    else
      name
  }

  # build the main tab content div
  contentDiv <- div(class=className("container"))
  if (!is.null(header))
    contentDiv <- tagAppendChild(contentDiv,div(class=className("row"), header))
  contentDiv <- tagAppendChild(contentDiv, tabset$content)
  if (!is.null(footer))
    contentDiv <- tagAppendChild(contentDiv,div(class=className("row"), footer))

  # build the page
  bootstrapPage(
    title = windowTitle,
    responsive = responsive,
    theme = theme,
    div(class=navbarClass,
      div(class="navbar-inner", containerDiv)
    ),
    contentDiv
  )
}
