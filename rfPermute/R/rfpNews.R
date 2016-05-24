#' @title Show the NEWS file
#' @description Show the \code{rfPermute} news file.
#' 
#' @export
#' 
rfpNews <- function () {
  file <- file.path(system.file(package = "rfPermute"), "NEWS.md")
  file.show(file)
}