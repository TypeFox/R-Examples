#' @name lazy.ref
#' @export lazy.ref
#' @export lazy.label
#' 
#' @title Reference Tables, Figures, Sections, and Pages
#' @description Provides the code to label and reference objects that 
#'   have a label statement with them
#'   
#' @param label A character(1) giving the name of the to be 
#'   created or referenced
#' @param text For HTML, the text to be hyperlinked for the reference.  
#'   If missing, this is set to \code{"(link)"}
#' @param page Indicates if the page number on which the label 
#'   lies should be returned or the object number.  This only applies 
#'   to LaTeX files
#' @param link for LaTeX files, should the reference link to the object
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' lazy.label("Label1")
#' lazy.ref("Label1")
#' lazy.ref("Label1", page=TRUE)
#' 

lazy.ref <- function(label, text, page=FALSE, link=TRUE){ 
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")

  #*** latex
  if (reportFormat == "latex"){
    code <- if (page)  paste("\\pageref{", label, "}", sep="") else paste("\\ref{", label, "}", sep="")
    if (link) code <- paste("\\hyperref[", label, "]{", code, "}", sep="")
  }
  
  if (reportFormat == "html"){
    if (missing(text)) text <- "(link)"
    code <- paste("<a href='#", label, "'> ", text, " </a>\n\n", sep="")
  }
  
  if (reportFormat == "markdown"){
    warning("labelling and referencing are not currently available in markdown")
    return("")
  }
  
  return(code)
}
