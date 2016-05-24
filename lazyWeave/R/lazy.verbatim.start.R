#' @name lazy.verbatim.start
#' @export lazy.verbatim.start
#' @export lazy.verbatim.end
#' 
#' @title Verbatim Environments
#' @description Text within a verbatim block appears exactly as type, 
#' including whitespace.  This is useful when inserting code into a document
#' 
#' @details A verbatim block takes any text entered and typsets it exactly as it was entered.  
#' White space is preserved and the font changes.  This is typically done to display
#' code, since the whitespace may preserve readability.
#' 
#' For HTML documents, this is done by opening a "<p ...>" tag with font family "monospace" and 
#' font "courier".  These are applicable until \code{lazy.verbatim.end} is called and the 
#' "</p>" tag is placed, closing the verbatim environment.
#' 
#' It should be noted that HTML code in this forced environment will still not render whitespace as in the
#' LaTeX verbatim environment.  This can be enforced by running the text in the environment through a function
#' like \code{gsub(" ", "&nbsp ", [text])} (\code{&nbsp} is the HTML character for a non-breaking space).
#' 
#' @author Benjamin Nutter
#' 

lazy.verbatim.start <- function(){
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  
  if (reportFormat == "latex") return("\\begin{verbatim}")

  #*** opens a paragraph, but does not end it.  Text inserted between lazy.verbatim.start and lazy.verbatim.end
  #*** will thus appear in monospace courier font
  if (reportFormat == "html") return("<p style='font-family:monospace, courier; font-size:11pt; text-align:left; font-style:none;
  font-weight:none; text-decoration:none;'>")
  
  if (reportFormat == "markdown") return("")
  
}
