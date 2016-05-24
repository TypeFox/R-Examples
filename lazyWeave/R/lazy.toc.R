#' @name lazy.toc
#' @export lazy.toc
#' 
#' @title Table of Contents and Other Lists
#' @description Designates the printing of a table of contents, list of 
#' figures, or list of tables.  Also provides functionality to manually edit
#' the contents of these lists
#' 
#' @param type Type of list to be printed or edited
#' @param add Determines if the list is printed or if an entry is added to the list
#' @param desc Only used when \code{add=TRUE}.  Gives the descriptive text of 
#'   the item being added to the list
#' @param withPage Determines if the page number of the entry is printed in the 
#'   table of contents. Only used when \code{add=TRUE}.
#' @param sec_unit Specifies the format for the new entry. For instance, will 
#' the new entry in the table of contents appear as a chapter, section, or 
#' subsection.  Used only when \code{withPage=TRUE}.
#' 
#' @details The level of detail a table of contents maintains is determined by 
#' the counter \code{tocdepth}.  In most cases, it is set to 3, giving chapter, 
#' section, and subsection.  To include subsubsections, it would be necessary to 
#' include \code{lazy.counter("tocdepth", value=4, fn="set")}.  
#' Use \code{value=5} to include paragraphs, and so forth.
#' 
#' @return Returns a string that designating that the table of contents is to 
#' be written, or an item to be added to a list.  This has no effect for HTML 
#' documents
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' lazy.toc()
#' lazy.toc("figures")
#' lazy.toc("tables", TRUE, "A brief description of the table")
#' lazy.toc("contents", TRUE, "Subsection 3", sec_unit="subsection")
#' 

lazy.toc <- function(type=c("contents", "figures", "tables"), add=FALSE, desc="",
    withPage=TRUE, sec_unit=c("chapter", "section", "subsection", "subsubsection", "part")){
  
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', 'markdown'")
  
  if (reportFormat == "latex"){

    fncall <- paste("%%", paste(deparse(match.call()), collapse=" "))
    type <- match.arg(type, c("contents", "figures", "tables"))
  
    sec_unit <- match.arg(sec_unit, c("chapter", "section", "subsection", "subsubsection", "part"))
  
    if (!add){
      code <- switch(type, 
          "contents" = "\\tableofcontents",
          "figures"  = "\\listoffigures",
          "tables"   = "\\listoftables")
    }
    else{
      code <- switch(type,
          "contents" = "toc",
          "figures" = "lof",
          "tables" = "lot")
      if (withPage) code <- paste("\\addcontentsline{", code, "}{", sec_unit, "}{", desc, "}", sep="")
      else code <- paste("\\addtocontents{", code, "}{", desc, "}", sep="")
    }
  
    code <- paste(fncall, "\n", code, "\n\n")
  }
  
  if (reportFormat %in% c("html", "markdown")){
    code <- ""
    warning("Tables of contents are not available in HTML or markdown reports.  Nothing has been done")
  }
  
  return(code)
}
