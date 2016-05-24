
#' Shiny interface for effectLite
#' 
#' This function calls a shiny interface for effectLite.
#' 
#' @param launch.browser Option will be passed on to \code{\link[shiny]{runApp}}
#' @export
effectLiteGUI <- function(launch.browser=TRUE){  
  shiny::runApp(system.file('elrshiny', package='EffectLiteR'),
                launch.browser=launch.browser)
}
