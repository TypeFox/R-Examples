#' Start lavaan.shiny
#' @title Launch lavaan.shiny User Interface
#' @return Nothing
#' @description lavaan.shiny() loads interactive user interface built using R shiny.
#' @details The interactive user interface is to provide an easy way for people who are learning how to work with the lavaan package and/or are not comfortable with the R command line system. Includes example data for testing out a few example analysis.
#' @keywords lavaan.shiny
#' @examples
#' \dontrun{
#' library(shiny)
#' lavaan.shiny()
#' }
#' @export
lavaan.shiny <- function() {

  shiny::runApp(appDir = system.file("shiny", package="lavaan.shiny"))

}
