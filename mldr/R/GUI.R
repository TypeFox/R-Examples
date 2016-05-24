#' The \pkg{mldr} package provides a basic, Shiny-based GUI to work with multilabel datasets.
#' You have to install the \pkg{shiny} package to be able to use this GUI.
#'
#' The user interface allows working with any of the previous loaded datasets, as well as loading
#' new ones. The GUI is structured into the following pages:
#' \itemize{
#'   \item{\strong{Main:}}{ This page is divided into two sections.
#'   The one at the left can be used to choose a previously loaded dataset,
#'   as well as to load datasets from files. The right part shows some basic
#'   statistics about the selected multilabel dataset.}
#'   \item{\strong{Labels:}}{ This page shows a table containing for each label its name, index, count, relative frequency
#'   and imbalance ratio (IRLbl). The page also includes a bar plot
#'   of the label frequency. The range of labels in the plot can be customized.}
#'   \item{\strong{Labelsets:}}{ This page shows a table containing for each labelset its representation and a counter.}
#'   \item{\strong{Attributes:}}{ This page shows a table containing for each attribute its name, type and a summary
#'   of its values.}
#'   \item{\strong{Concurrence:}}{ This page shows for each label the number of instances in which it appears and its
#'   mean SCUMBLE measure, along with a plot that shows the level of concurrence among the selected labels. Clicking the
#'   labels in the table makes it possible to add/remove them from the plot.}
#' }
#'
#' The tables shown in these pages can be sorted by any of its fields, as well as filtered. The content of the tables
#' can be copied to clipboard, printed and saved in CSV and Microsoft Excel format.
#'
#' @title Launches the web-based GUI for mldr
#' @return Nothing
#' @description Loads an interactive user interface in the web browser, built using R shiny.
#' @examples
#'\dontrun{
#' library(mldr)
#' mldrGUI()
#'}
#' @import shiny
#' @export

mldrGUI <- function() {
  shiny::runApp(appDir = system.file("shiny", package="mldr"), launch.browser = TRUE)

  invisible()
}
