#' web.embryogrowth runs a web application for basic functions of embryogrowth
#' @title Run a web application for basic functions of embryogrowth
#' @author Marc Girondot
#' @author Maria Sousa Martins
#' @return Nothing
#' @description Run a shiny application for basic functions of embryogrowth
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' .web.embryogrowth()
#' }
#' @export


.web.embryogrowth <- function() {

  
  if (requireNamespace("shiny", quietly = TRUE)) {
    shiny::runApp(appDir = system.file("shiny", package="embryogrowth"),launch.browser =TRUE)
  } else {
    warning("shiny package is required for this function")
  }
  

}
