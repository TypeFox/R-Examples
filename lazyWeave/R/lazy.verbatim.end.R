#' @describeIn lazy.verbatim.start
#' 

lazy.verbatim.end <- function(){
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  if (reportFormat == "latex") return("\\end{verbatim}")
  
  if (reportFormat == "html") return("</p>")
  
  if (reportFormat == "markdown") return("")
}
