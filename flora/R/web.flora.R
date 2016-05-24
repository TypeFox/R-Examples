#' Web front end
#' 
#' This function starts a local webserver to run the shiny app distributed with
#' the package.
#' 
#' @export
web.flora <- function() {
  message("Press escape at any time to stop the application.\n")
  runApp(system.file("plantminer", package = "flora"))
}