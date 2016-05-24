#' Run npregfast example
#'
#' Launch a Shiny app that shows a demo of what can be done with
#' the package.
#'
#' This example is also
#' \href{http://sestelo.shinyapps.io/npregfast/}{available online}.
#'
#' @examples
#' ## Only run this example in interactive R sessions
#' if (interactive()) {
#'   runExample()
#' }
#' @importFrom shinyjs colourInput useShinyjs
#' @importFrom wesanderson wes_palettes
#' @export



#' @export

runExample <- function() {
 
  if (!requireNamespace("npregfast", quietly = TRUE)) {
    stop('`npregfast` package is required for this function.\nPlease install it with `install.packages("npregfast")`',
         call. = FALSE)
  }
  
  appDir <- system.file("shiny_examples", "demo", package = "npregfast")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `npregfast`.",
         call. = FALSE)
  }
  
  
  shiny::runApp(appDir, display.mode = "normal")
  
}
