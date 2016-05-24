#' Shiny Visualisation App
#'
#' Provides an interactive web app for the running the IncucyteDRC workflow
#'
#' @return Launches an interactive Shiny application
#' @export
#' @import shiny
shinyVisApp <- function() {

    shiny::shinyApp(
        ui = shinyVisUI(),
        server = function(input, output) {
            shinyVisServer(input, output)
        }
    )
}
