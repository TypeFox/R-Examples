#' @name lazy.footnote
#' @export lazy.footnote
#' 
#' @title Add a Footnote
#' @description Adds a footnote to the text in a .tex file.  
#' In HTML files, an endnote is produced.  Links are 
#' established between the text and the endnote for user convenience
#' 
#' @param text Text for the footnote
#' @param number Footnote number
#' @param translate Determines if \code{\link{latexTranslate}}
#'  is applied to \code{text}.
#' @param name For HTML, a reference name to the endnote
#' @param ref For HTML, a reference name to go back to the text (from the endnote).
#' @param counter For HTML, counter to use for numbering the endnotes
#' @param size For HTML, the text size to use for the endnote
#' 
#' @details Be warned that this is not a perfect function on the LaTeX side.
#'   If you use special characters that require that \code{latexTranslate} be 
#'   turned off, you'll also need to turn of \code{latexTranslate} in
#'    \code{lazy.write}.  I'm including this as is for ease of use, but it 
#'    does have some kinks to work out.
#'    
#'    \code{name} and \code{ref} are used to create links between the 
#'    footnote marking in the text and the actual footnote.  Clicking
#'    on the links created will allow the reader to go back and forth 
#'    between them.  The names may be similar, but exact matches may
#'    confuse the browser.
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' lazy.footnote("Enter a reference to an article in this argument", number=3)
#' lazy.footnote(lazy.citation(), number=4)

lazy.footnote <- function(text, number=NULL, translate=FALSE,
                          name, ref, counter="footnote", size=8){
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  if (reportFormat == "latex"){
    num <- if (is.null(number)) "" else paste("[", number, "]", sep="")  
    code <- paste("\\footnote", num, "{", text, "}", sep="")
  }
  
  if (reportFormat == "html"){
    if (!is.null(number)){
      if (!is.numeric(number)) stop("'number' must be numeric")
      lazy.counter(counter, number, fn="set")
    } 
    
    end.value <- lazy.counter(counter, fn="value")
    lazy.counter(counter, end.value + 1, "set")
    
    code <- paste("<sup>[<a name='", name, "' href='#", ref, "'>", end.value, "</a>]</sup>", sep="")
    
    to.add <- lazy.text(paste("<sup>[<a name='", ref, "' href='#", name, "'>", end.value, "</a>]</sup>", sep=""), text, size=size)
    assign("HTML.FOOTNOTES", paste(get("HTML.FOOTNOTES", envir=options()$htmlCounters), to.add, sep="\n"), envir=options()$htmlCounters) 
  }
  
  if (reportFormat == "markdown"){
    if (!is.null(number)){
      if (!is.numeric(number)) stop("'number' must be numeric")
      lazy.counter(counter, number, fn="set")
    } 
    
    end.value <- lazy.counter(counter, fn="value")
    lazy.counter(counter, end.value + 1, "set")
    
    code <- end.value
    to.add <- paste("^", end.value, "^ ", text, sep="")
    assign("HTML.FOOTNOTES", paste(get("HTML.FOOTNOTES", envir=options()$htmlCounters), to.add, sep="\n\n"), envir=options()$htmlCounters)
  }   
  
  return(code)
}
