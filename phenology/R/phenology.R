#' phenology runs a shiny application for basic functions of phenology
#' @title Run a shiny application for basic functions of phenology
#' @author Marc Girondot
#' @return Nothing
#' @description Run a shiny application for basic functions of phenology
#' @examples
#' \dontrun{
#' library(phenology)
#' phenology()
#' }
#' @export


phenology <- function() {

runApp(appDir = system.file("shiny", package="phenology"),launch.browser =TRUE)

}
