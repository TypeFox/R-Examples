
#' @title Run statnetWeb
#'
#' @description Runs the statnetWeb shiny application, a GUI for the statnet
#' suite of network analysis packages.
#'
#' @keywords GUI, shiny
#' @export
#' @examples
#' \dontrun{
#' run_sw()
#' }
#'
#'
run_sw <- function() {
  shiny::runApp(system.file("shiny", package = "statnetWeb"))
}