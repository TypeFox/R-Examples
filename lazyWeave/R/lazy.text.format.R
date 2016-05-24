#' @name lazy.text.format
#' @export lazy.text.format
#' 
#' @title Format Text
#' @description Applies italic, bold, or underlining to a piece of text.  
#'   May be used within \code{lazy.text} to add emphasis when the entire 
#'   paragraph does not need to be formatted
#' 
#' @param text Text to be formatted
#' @param italic Logical.  Specifies if text should be italic
#' @param bold Logical.  Specifies if text should be bold
#' @param underline Logical. Specifies if text should be underlined
#' @param translate Logical. Specifies if text should be passed through 
#'   \code{latexTranslate} before printing
#'   
#' @details This function differs from \code{lazy.text} in that 
#' \code{lazy.text} produces a paragraph of formatted text while
#' \code{lazy.text.format} produces smaller blocks.  This allows for 
#' smaller bits of text to be formatted for emphasis
#' (see the last example below).
#' 
#' @author Benjamin Nutter
#' 
#' @examples 
#' lazy.text.format("This is the text")
#' lazy.text.format("This is the text", italic=TRUE)
#' lazy.text.format("This is the text", bold=TRUE)
#' lazy.text.format("This is the text", italic=TRUE, bold=TRUE)
#' 
#' lazy.text("The percentage of defective lightbulbs in this sample was ", 
#'           lazy.text.format("30\%", italic=TRUE),
#'           ". Clearly, this is unacceptable.")
#' 


lazy.text.format <- function(text, italic=FALSE, bold=FALSE, 
                             underline=FALSE, translate=TRUE){

  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  if (reportFormat == "latex"){
    if (translate) text <- latexTranslate(text)
  
    if (underline) text <- paste("\\ul{", text, "}", sep="")
    if (bold)      text <- paste("\\textbf{", text, "}", sep="")
    if (italic)    text <- paste("\\emph{", text, "}", sep="")
  }
  
  if (reportFormat == "html"){
    if (underline) text <- paste("<ul>", text, "</ul>")
    if (italic) text <- paste("<i>", text, "</i>")
    if (bold) text <- paste("<b>", text, "</b>")  
  }
  
  if (reportFormat == "markdown"){
    if (italic) text <- paste0("_", text, "_")
    if (bold) text <- paste0("**", text, "**")
  }
  
  return(text)
}
