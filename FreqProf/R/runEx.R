#' Run interactive FreqProf example (Shiny App)
#'
#' @export
#' @examples
#' \donttest{
#'  runEx()
#' }
runEx <- function() {
  appDir <- system.file("shinyapp", package = "FreqProf")
  if(appDir == "") {
    stop("Could not find example directory. Try reinstalling `FreqProf`.",
         call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}