#' @title Examples from the SciencesPo Package
#'
#' @description Launchs a Shiny app that shows a demo.
#' @param example The name of the shiny application.
#' @examples
#' # A demo of what \code{\link[SciencesPo]{politicalDiversity}} does.
#' if (interactive()) {
#' runExample(politicalDiversity)
#' }
#' @export
#'
`runExample` <- function(example){
# locate all the shiny app examples that exist
validExamples <- list.files(system.file("shiny-examples", package = "SciencesPo"))
validExamplesMsg <-
paste0(
"Valid examples are: '",
paste(validExamples, collapse = "', '"),
"'")
# if an invalid example is given, throw an error
 if (missing(example) || !nzchar(example) ||
   !example %in% validExamples) {
   stop(
'Please run `runExample()` with a valid example app as an argument.\n',
validExamplesMsg,
     call. = FALSE)
}
# find and launch the app
appDir <- system.file("shiny-examples", example, package = "SciencesPo")
shiny::runApp(appDir, display.mode = "normal")
}##--end of runExample
NULL
