#' Start IRTShiny
#' @title This function will start IRTShiny
#' @return Nothing
#' @description An interactive Shiny application for running a IRT analysis.
#' @details This starts the IRT Shiny application on the users local computer. 
#' @keywords IRT
#' @examples
#' \dontrun{
#' library(shiny)
#' library(shinyAce)
#' library(psych)
#' library(CTT)
#' library(ltm)
#' library(beeswarm)
#' library(parallel)
#' startIRT()
#' }
#' @export

startIRT <- function() {
  
  shiny::runApp(appDir = system.file("IRT", package="IRTShiny"))
  
}
