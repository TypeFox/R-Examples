#' @name lazy.page.break
#' @export lazy.page.break
#' 
#' @title Start New Page in LaTeX
#' @description Insert code to start a new page in a document
#' 
#' @details For HTML documents, page breaks will not show up in the browser, 
#' but will show up when printed. 
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' \dontrun{
#' lazy.write(
#' lazy.file.start(),
#' lazy.text("First we type something on the first page"),
#' lazy.page.break(),
#' lazy.text("Then we type something on the second page"),
#' lazy.file.end(),
#' OutFile="Example 1.tex")
#' 
#' unlink("Example 1.tex")
#' }


lazy.page.break <- function(){
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  if (reportFormat == "latex") return("%% lazy.page.break()\n\\newpage")
  if (reportFormat == "html") return("<!--html_page_break()-->\n<div style='page-break-before:always'></div>")
  if (reportFormat == "markdown") return("*******")
}
