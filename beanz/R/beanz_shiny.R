
##-----------------------------------------------------------------------------------
##                     run Shiny Gui
##-----------------------------------------------------------------------------------
#' Run Web-Based BEANZ application
#'
#' Call Shiny to run \code{beanz} as a web-based application
#'
#'
#'
#' @export
#'
run.beanz <- function() {
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Shiny needed for this function to work. Please install it.",
             call. = FALSE)
    }

    appDir <- system.file("shiny", package = "beanz")
    if (appDir == "") {
        stop("Could not find Shiny directory. Try re-installing `beanz`.",
             call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal");
}
