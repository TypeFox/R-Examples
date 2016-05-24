#' Create a navigation list panel
#'
#' Create a navigation list panel that provides a list of links on the left
#' which navigate to a set of tabPanels displayed to the right.
#'
#' @param ... \code{\link[shiny]{tabPanel}} elements to include in the navlist
#' @param id If provided, you can use \code{input$}\emph{\code{id}} in your
#'   server logic to determine which of the current navlist items is active. The
#'   value will correspond to the \code{value} argument that is passed to
#'   \code{\link[shiny]{tabPanel}}.
#' @param selected The \code{value} (or, if none was supplied, the \code{title})
#'   of the navigation item that should be selected by default. If \code{NULL},
#'   the first navigation will be selected.
#' @param well \code{TRUE} to place a well (gray rounded rectangle) around the
#'   navigation list.
#' @param fluid \code{TRUE} to use fluid layout; \code{FALSE} to use fixed
#'   layout.
#' @param widths Column withs of the navigation list and tabset content areas
#'   respectively.
#'
#' @details You can include headers within the \code{navlistPanel} by
#' including plain text elements in the list; you can include separators by
#' including "------" (any number of dashes works).
#'
#' @examples
#' library(shiny)
#'
#' fluidPage(
#'
#'   titlePanel("Application Title"),
#'
#'   navlistPanel(
#'     "Header",
#'     tabPanel("First"),
#'     tabPanel("Second"),
#'     "-----",
#'     tabPanel("Third")
#'   )
#' )
#' @export
navlistPanel <- function(...,
                         id = NULL,
                         selected = NULL,
                         well = TRUE,
                         fluid = TRUE,
                         widths = c(4, 8)) {

  # text filter for defining separtors and headers
  textFilter <- function(text) {
    if (grepl("^\\-+$", text))
      tags$li(class="divider")
    else
      tags$li(class="nav-header", text)
  }

  # build the tabset
  tabs <- list(...)
  tabset <- buildTabset(tabs,
                        "nav nav-list",
                        textFilter,
                        id,
                        selected)

  # create the columns
  columns <- list(
    column(widths[[1]], class=ifelse(well, "well", ""), tabset$navList),
    column(widths[[2]], tabset$content)
  )

  # return the row
  if (fluid)
    fluidRow(columns)
  else
    fixedRow(columns)
}


buildTabset <- function(tabs,
                        ulClass,
                        textFilter = NULL,
                        id = NULL,
                        selected = NULL) {

  # build tab nav list and tab content div

  # add tab input sentinel class if we have an id
  if (!is.null(id))
    ulClass <- paste(ulClass, "shiny-tab-input")

  tabNavList <- tags$ul(class = ulClass, id = id)
  tabContent <- tags$div(class = "tab-content")
  firstTab <- TRUE
  tabsetId <- shiny:::p_randomInt(1000, 10000)
  tabId <- 1
  for (divTag in tabs) {

    # check for text; pass it to the textFilter or skip it if there is none
    if (is.character(divTag)) {
      if (!is.null(textFilter))
        tabNavList <- tagAppendChild(tabNavList, textFilter(divTag))
      next
    }

    # compute id and assign it to the div
    thisId <- paste("tab", tabsetId, tabId, sep="-")
    divTag$attribs$id <- thisId
    tabId <- tabId + 1

    tabValue <- divTag$attribs$`data-value`
    if (!is.null(tabValue) && is.null(id)) {
      stop("tabsetPanel doesn't have an id assigned, but one of its tabPanels ",
           "has a value. The value won't be sent without an id.")
    }


    # function to append an optional icon to an aTag
    appendIcon <- function(aTag, iconClass) {
      if (!is.null(iconClass)) {
        # for font-awesome we specify fixed-width
        if (grepl("fa-", iconClass, fixed = TRUE))
          iconClass <- paste(iconClass, "fa-fw")
        aTag <- tagAppendChild(aTag, shiny::icon(name = NULL, class = iconClass))
      }
      aTag
    }

    # check for a navbarMenu and handle appropriately
    if (inherits(divTag, "shiny.navbarmenu")) {

      # create the a tag
      aTag <- tags$a(href="#",
                     class="dropdown-toggle",
                     `data-toggle`="dropdown")

      # add optional icon
      aTag <- appendIcon(aTag, divTag$iconClass)

      # add the title and caret
      aTag <- tagAppendChild(aTag, divTag$title)
      aTag <- tagAppendChild(aTag, tags$b(class="caret"))

      # build the dropdown list element
      liTag <- tags$li(class = "dropdown", aTag)

      # build the child tabset
      tabset <- buildTabset(divTag$tabs, "dropdown-menu")
      liTag <- tagAppendChild(liTag, tabset$navList)

      # don't add a standard tab content div, rather add the list of tab
      # content divs that are contained within the tabset
      divTag <- NULL
      tabContent <- tagAppendChildren(tabContent,
                                      list = tabset$content$children)
    }
    # else it's a standard navbar item
    else {
      # create the a tag
      aTag <- tags$a(href=paste("#", thisId, sep=""),
                     `data-toggle` = "tab",
                     `data-value` = tabValue)

      # append optional icon
      aTag <- appendIcon(aTag, divTag$attribs$`data-icon-class`)

      # add the title
      aTag <- tagAppendChild(aTag, divTag$attribs$title)

      # create the li tag
      liTag <- tags$li(aTag)
    }

    if (is.null(tabValue)) {
      tabValue <- divTag$attribs$title
    }

    # If appropriate, make this the selected tab (don't ever do initial
    # selection of tabs that are within a navbarMenu)
    if ((ulClass != "dropdown-menu") &&
       ((firstTab && is.null(selected)) ||
        (!is.null(selected) && identical(selected, tabValue)))) {
      liTag$attribs$class <- "active"
      divTag$attribs$class <- "tab-pane active"
      firstTab = FALSE
    }

    divTag$attribs$title <- NULL

    # append the elements to our lists
    tabNavList <- tagAppendChild(tabNavList, liTag)
    tabContent <- tagAppendChild(tabContent, divTag)
  }

  list(navList = tabNavList, content = tabContent)
}
