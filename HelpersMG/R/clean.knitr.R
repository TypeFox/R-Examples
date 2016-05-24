#' clean.knitr deletes temporary files created during knitr compile
#' @title Delete temporary files created during knitr compile
#' @author Marc Girondot
#' @return Nothing
#' @description Delete temporary files created during knitr compile in 
#' working directory.\cr
#' This function works only in UNIX system (LINUX or MacOSX).\cr
#' @examples
#' \dontrun{
#' clean.knitr()
#' }
#' @export


clean.knitr <- function() {
  if (.Platform$OS.type=="unix") {
    
  system(paste0("cd '", getwd(), "';rm *.gz;rm *.toc;rm *.tex;rm *.log"))
  } else {
    warning("This function works only in UNIX systems")
  }
}
