#' @name setHtmlOptions
#' @export setHtmlOptions
#' 
#' @title lazyWeave HTML Report Options
#' @description Sets or changes options for report counters, font, font size, 
#'   and footnotes.  The counter options apply to HTML and Markdown documents,
#'   and the font options only apply to HTML.
#' 
#' @param table Sets the HTML table counter
#' @param figure Sets the HTML figure counter
#' @param footnote Sets the HTML footnote counter
#' @param chapter Sets the HTML chapter counter
#' @param section Sets the HTML section counter
#' @param subsection Sets the HTML section counter
#' @param subsubsection Sets the HTML subsubsection counter
#' @param font.family Sets the HTML font family
#' @param font Sets the HTML font
#' @param font.size Sets the HTML font size
#' 
#' @details
#' For all arguments, a value of \code{NULL} results in no action.
#'
#' The HTML counters are used to maintain a somewhat consistent appearance 
#' between HTML and LaTeX reports.  Since HTML doesn't have counters, a series 
#' of variables is inserted into the Global Environment.  These are hidden 
#' from view and are incremented automatically by   \code{lazyWeave} function.  
#' Manipulating these variables directly is strongly discouraged.  They can 
#' all be managed by \code{lazy.counter}.
#' 
#' \code{setHtmlOptions} is a convenience function that can change all of the 
#' global variables simultaneously (as opposed to repeated calls to 
#' \code{lazy.counter}).  However, this is the recommended way to change 
#' font family, font, and font size.
#' 
#' To change the report format, use the code 
#' \code{options(lazyReportFormat = "latex")}, 
#' \code{options(lazyReportFormat = "html")} or,
#' \code{options(lazyReportFormat = "markdown")}
#' 
#' @author Benjamin Nutter
#' 

setHtmlOptions <- function(table=NULL, figure=NULL, footnote=NULL,
                         chapter=NULL, section=NULL, subsection=NULL,
                         subsubsection=NULL, 
                         font.family=NULL, font=NULL, font.size=NULL){
  if (!is.null(table)) assign("HTML.COUNTER.TABLE", table, envir=options()$htmlCounters)
  if (!is.null(figure)) assign("HTML.COUNTER.FIGURE", figure, envir=options()$htmlCounters)
  if (!is.null(footnote)) assign("HTML.COUNTER.FOOTNOTE", footnote, envir=options()$htmlCounters)
  if (!is.null(chapter)) assign("HTML.COUNTER.CHAPTER", chapter, envir=options()$htmlCounters)
  if (!is.null(section)) assign("HTML.COUNTER.SECTION", section, envir=options()$htmlCounters)
  if (!is.null(subsection)) assign("HTML.COUNTER.SUBSECTION", subsection, envir=options()$htmlCounters)
  if (!is.null(subsubsection)) assign("HTML.COUNTER.SUBSUBSECTION", subsubsection, envir=options()$htmlCounters)
  if (!is.null(font.family)) assign("HTML.FONT.FAMILY", font.family, envir=options()$htmlCounters)
  if (!is.null(font)) assign("HTML.FONT.FONT", font, envir=options()$htmlCounters)
  if (!is.null(font.size)) assign("HTML.FONT.SIZE", font.size, envir=options()$htmlCounters)
}
