#' @name lazy.landscape
#' @export lazy.landscape
#' 
#' @title Landscape Page Orientation
#' @description Allows the user to turn on or off landscape page orientation
#' 
#' @param begin If \code{TRUE}, begins the landscape environment.  
#'   Otherwise, the environment is closed
#'   
#' @details This function has no effect on HTML or markdown files.  
#'   The page orientation is determined by the browser
#' 
#' @return Returns a string that either begins or ends a 
#'  landscape environment
#'  
#' @author Benjamin Nutter
#' 
#' @examples 
#' lazy.landscape()
#' lazy.landscape(begin=FALSE)

lazy.landscape <- function(begin=TRUE){
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")

  if (reportFormat == "latex") if (begin) return("\\begin{landscape}\n") else return("\\end{landscape}")

  if (reportFormat == "html") return("")
  
  if (reportFormat == "markdown"){
    warning("Landscape orientation is not available for markdown")
    return("")
  }
}
