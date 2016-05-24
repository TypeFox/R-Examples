#' @export
runShinyPIPE <- function() {
  appDir <- system.file("shiny", "shinyPIPE", package = "pipe.design")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `pipe.design`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}