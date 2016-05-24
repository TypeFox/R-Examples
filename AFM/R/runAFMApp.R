#' Launch the AFM shiny application
#'
#' Launch the AFM shiny graphical interface to access most of the fonctionalities of the AFM library
#'
#' @export
#' @author M.Beauvais
#' @examples
#' \dontrun{
#' install.packages("AFM")
#' AFM::runAFMapp()
#' }
runAFMApp <- function() {
  appDir <- system.file("shiny", "AFM-desktop", package = "AFM")
  if (appDir == "") {
    stop("Could not find directory. Try re-installing `AFM`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}