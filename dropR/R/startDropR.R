#' Start the DropR Shiny App
#' 
#' Starts the interactive web application to use dropR in your web browser. 
#' Make sure to use Google Chrome or Firefox for best experience.
#' 
#' @examples
#' \dontrun{startdropR()}
#' @export
startDropR <- function(){
  shiny::runApp(system.file('dropR_shiny', package='dropR'),launch.browser = T)  
}
