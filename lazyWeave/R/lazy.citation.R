#' @name lazy.citation
#' @export lazy.citation
#' 
#' @title Add R or R Package Citation
#' @description Generates code for the citation of R or an R package.
#' 
#' @param pkg a character(1) vector giving the name of a package.  If \code{NULL}, a citation for R is produced.
#' @param author Include author name
#' @param title Include title of package
#' @param org Include organization name
#' @param address Include address
#' @param volume Include volume
#' @param year include year of publication
#' @param note include the note on the citation.
#' 
#' @details Not every option is populated in every package.  
#' Future improvements might include automatic detection of NULL fields, 
#' but for now, observing the output with all the options set to \code{TRUE}
#' will tell you which ones are empty.
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' lazy.citation()
#' lazy.citation(pkg="lazyWeave", org=FALSE, address=FALSE, volume=FALSE)

lazy.citation <- function(pkg=NULL, author=TRUE, title=TRUE, org=TRUE, 
                          address=TRUE, volume=TRUE, year=TRUE, note=TRUE){
  
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  

  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")
  

  #*** get the right left quote characters for the report format
  quote.string <- if (reportFormat == "latex") "``"              
  else if (reportFormat %in% c("html", "markdown")) "\""
  
  #*** Construct the citation
  cit <- if (is.null(pkg)) utils::citation() else utils::citation(pkg)

  paste( if (author) paste(paste(cit$author, collapse=", "), ", ", sep="") else "",
         if (title)  paste(quote.string, cit$title, ",\" ", sep="") else "",
         if (org)    paste(cit$organization, ", ", sep="") else "",
         if (address) paste(cit$address, ", ", sep="") else "",
         if (volume)  paste("Vol. ", cit$volume, ", ", sep="") else "",
         if (year)    paste("(", cit$year, ") ", sep="") else "",
         if (note)    paste(cit$note, ".", sep="") else "", sep="")
}

