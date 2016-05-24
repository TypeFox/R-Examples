#' Launch ggraptR in the default browser
#'
#' @details See \url{http://github.com/cargomoose/raptR} for documentation and tutorials
#'
#' @examples
#' if (interactive()) {
#'   ggraptR()
#' }
#' @import ggplot2
#' @import DT
#' @importFrom dplyr group_by_
#' @importFrom dplyr summarise
#' @importFrom dplyr summarise_each
#' @import futile.logger
#' @import ggthemes
#' @import shinyBS
#' @import shinyjs
#' @export
ggraptR <- function() {
  appDir <- system.file("ggraptR", package = "ggraptR")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}

