#' Start weightr in Shiny
#'
#' This function allows you to launch the Shiny application locally.
#' @keywords weightr
#' @export
#' @importFrom stats model.matrix optim pchisq pnorm qnorm
#' @examples
#' \dontrun{
#' library(shiny)
#' library(shinyBS)
#' shiny_weightr()
#' }

shiny_weightr <- function() {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The R package 'shiny' is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("shinyBS", quietly = TRUE)) {
    stop("The R package 'shinyBS' is needed for the application to display correctly. Please install it.",
         call. = FALSE)
  }

  shiny::runApp(appDir = system.file("shiny", package="weightr"))

}
