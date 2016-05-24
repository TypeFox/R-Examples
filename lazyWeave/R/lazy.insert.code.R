#' @name lazy.insert.code
#' @export lazy.insert.code
#' 
#' @title Insert Code to a TeX or HTML document
#' @description Places TeX code from a file into a TeX document.  
#' Text is placed within a verbatim environment.  For HTML documents, 
#' a verbatim-like environment is created
#' 
#' @param file Filename containing the code to be inserted
#' @param prompt.symbol Character to be placed at the left most side 
#'   of the page.  Defaults to the system's current prompt symbol
#' @param lines A vector giving the lines in \code{file} to be 
#'   inserted into the document
#'   
#' @details Text is inserted in a verbatim environment to preserve 
#'   whitespace.  This function is performed better by \code{Sweave} 
#'   and \code{knitr}.  Those packages will also display any result 
#'   printed by the code.  This function will not display results.
#'   
#'   With HTML, the font family is set to monospace, and the font is 
#'   set to courier.  All of the spaces are replaced with a 
#'   non-breaking space to preserve the white space in the code.
#'   
#' @author Benjamin Nutter
#' 

lazy.insert.code <- function(file, prompt.symbol=options()$prompt, lines){
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  collapse.char = if (reportFormat == "latex") "\n" else if (reportFormat == "html") "<br>" else "\n\n"
  
  file <- readLines(file)
  if (missing(lines)) lines <- 1:length(file)
  if (reportFormat=='markdown') file <- paste("`", file, "`", sep="")
  file <- paste(paste(prompt.symbol, file[lines], collapse=collapse.char))
  if (reportFormat == "html") file <- gsub(" ", "&nbsp ", file)
  paste(lazy.verbatim.start(), file, lazy.verbatim.end(), sep="\n")
}
