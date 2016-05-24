#' Graphical user interface to \pkg{RJafroc} functions
#' 
#' Start a graphical user interface which is functionally similar to Windows JAFROC software (Version 4.2.1)
#' 
#' @param useBrowser a logical variable, if TRUE the default internet browser is used, otherwise RStudio's
#' internal browser is used. Default is \code{FALSE}. See "Details".
#' 
#' @details For Windows users, we suggest setting \code{useBrowser} to \code{TRUE} due to a bug in RStudio's (Version 0.99.467) 
#' internal broswer which prevents saving plots.
#' 
#' @examples
#' \dontrun{
#' ## For Windows users:
#' RJafrocGui(useBrowser = TRUE)
#' 
#' ## For other users:
#' RJafrocGui()
#' }
#' 
#' @export
#' 
#' @import shiny
RJafrocGui <- function(useBrowser = FALSE){  
  appDir <- system.file("GUI", package = "RJafroc")
  if (useBrowser){
    runApp(appDir, launch.browser = useBrowser)
  }else{
    runApp(appDir)
  }
}