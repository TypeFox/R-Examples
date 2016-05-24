#' @name lazy.file.end
#' @export lazy.file.end
#' 
#' @title End LaTeX Documents
#' @description Provides the code to end a LaTeX document
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' lazy.file.end()
#' 
lazy.file.end <- function(){

  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  if (reportFormat == "latex") return("\n\n\\end{document}")
  else if (reportFormat == "html") return("\n\n</html>\n")
  else if (reportFormat == "markdown") return("")
}

