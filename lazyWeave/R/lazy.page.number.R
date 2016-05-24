#' @name lazy.page.number
#' @export lazy.page.number
#' 
#' @title Change Page Numbering
#' @description Allows page numbering style to be changed.  For instance, 
#' from roman numerals in an introduction to arabic numbers in the body of a 
#' report
#' 
#' @param num_style A character(1) giving the numbering style for the page
#' 
#' @details Each time the page numbering is changed, the page counter is 
#' reset to 0 (meaning the next page will be numbered 1).  If the number
#' needs to be preserved, this can be done using \code{lazy.counter} with
#'  \code{counter="page"}.
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' lazy.page.number("Roman")
#' 

lazy.page.number <- function(num_style = c("arabic", "roman", "Roman", "alph", "Alph")){
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', 'markdown'")
  
  if (reportFormat == "latex"){
    num_style <- match.arg(num_style, c("arabic", "roman", "Roman", "alph", "Alph"))
    return(paste("\\pagenumbering{", num_style, "}\n\n", sep=""))
  }
  if (reportFormat == "html") return("")
  
  if (reportFormat == "markdown") return("")
    
}
