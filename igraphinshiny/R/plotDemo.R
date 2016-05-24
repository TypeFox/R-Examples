#' @export

plotDemo <- function() {
  appDir <- system.file("igraphPlotDemo", package = "igraphinshiny")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `igraphinshiny`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
