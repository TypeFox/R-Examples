#' @name lazy.link
#' @export lazy.link
#' 
#' @title Links to Webpages or External Documents
#' @description While \code{lazy.ref} provides the option of linking 
#' to areas within a document, \code{lazy.link} provides the option of 
#' linking to areas outside of the document.  Web pages are perhaps 
#' the most obvious example, but links could also go to files on a 
#' directory
#' 
#' @param url A character(1) giving the URL for the link or a file path
#' @param text The text to be highlighted as the link.  If this 
#'   is missing, \code{url} is used
#' @param web When \code{TRUE}, \code{"http://"} is added to 
#'   \code{url}, (if not already present), to ensure a link to the 
#'     web.  For files on a local dis, set this to \code{FALSE}
#' @param secure Should the link be to a secure site "https://".
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' lazy.link("https://github.com/nutterb/lazyWeave", secure=TRUE)
#' 

lazy.link <- function(url, text, web=TRUE, secure=FALSE){
  
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")
  
  #*** append http (or https) if not on url
  if (web){
    first.seven <- substr(url, 1, if(secure) 8 else 7)
    if (!first.seven %in% if (secure) "https://" else "http://") 
      url <- paste(if (secure) "https://" else "http://", url, sep="")
  }
  
  if (reportFormat == "latex"){
    code <- if (missing(text)) paste("\\url{", url, "}\n", sep="") else paste("\\href{", url, "}{", text, "}\n", sep="")
  }
  
  if (reportFormat == "html"){
    code <- paste("<a href='", url, "'>", if (missing(text)) url else text, " </a>\n\n", sep="")
  }
  
  if (reportFormat == "markdown"){
    fncall <- ""
    code <- paste("[", if (missing(text)) url else text, "](", url, ")", sep="")
  }
 
  return(paste(fncall, code))
}
