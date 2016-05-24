#' @name lazy.text
#' @export lazy.text
#' 
#' @title Paragraphs in LaTeX
#' @details Write paragraphs in LaTeX code
#' 
#' @param ... Text and other objects to be included in the paragraph
#' @param title A title for the paragraph
#' @param align alignment of the paragraph.  Options are \code{"left", 
#'   "center", "right"}.
#' @param italic Logical.  Indicates if the text should be italicized
#' @param bold Logical.  Indicates if the text shoudl be bolded
#' @param underline Logical.  Indicates if the text should be underlined
#' @param sep Character.  Denotes the separation string for the items in 
#'   \code{...} when they are pasted together
#' @param translate Toggles if inputs in \code{x} should be passed through 
#'   \code{latexTranslate}.  This should be set to \code{FALSE} if writing
#'   custom code.
#' @param font HTML font for the paragraph. Defaults to the HTML option 
#'   (see \code{\link{setHtmlOptions}}).
#' @param family HTML font family for the paragraph. Defaults to the HTML 
#'   option (see \code{\link{setHtmlOptions}}).
#' @param size Text size of the paragraph.   Defaults to the HTML option 
#'   (see \code{\link{setHtmlOptions}}). May be an integer or a LaTeX size 
#'   descriptor. See "Details" for options
#'   
#' @details Options for text size are
#' \tabular{ll}{
#'   tiny \tab smallest \cr
#'   scriptsize \tab \cr
#'   footnotesize \tab \cr
#'   small \tab \cr
#'   normalsize \tab \cr
#'   large \tab \cr
#'   Large \tab \cr
#'   LARGE \tab \cr
#'   huge \tab \cr
#'   Huge \tab BIGGEST\cr
#' }
#' 
#' When size is denoted as an integer, as necessary for HTML, it is mapped to a LaTeX size using 
#' \code{map.size}.  Likewise, the LaTeX descriptors can be mapped to integers using \code{map.size}.
#' 
#' Additional formatting may be applied using commands such as \code{textbf\{\}} for bold text, \code{emph\{\}} for italics, and
#' \code{ul\{\}} for underlined text (assuming the \code{soul} package is available), but doing so is probably
#' easier using \code{lazy.text.format}.   
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' \dontrun{
#' lazy.write(
#'   lazy.file.start(),
#'   lazy.text("Typically we want our paragraphs to be left 
#'     justified.  This is often what we expect to see when reading."),
#'   lazy.text("However, we may also have occasions where we would 
#'     like to center our text.  It's one of many ways we can make the 
#'     words stand out on the page", align="center"),
#'   lazy.text("A more traditional way to make the text stand out might be
#'     to use bold text or italics, such as these", bold=TRUE, italic=TRUE),
#'   lazy.file.end(),
#'   OutFile="Example 1.tex")
#'   
#' unlink("Example 1.tex")
#' }
#' 

lazy.text <-
function(..., title=NULL, align="left",
    italic=FALSE, bold=FALSE, underline=FALSE, sep="",
    translate=TRUE, font, family, size){
    
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")
  
  if (missing(font)) font <- get("HTML.FONT.FONT", envir=options()$htmlCounters)
  if (missing(family)) family <- get("HTML.FONT.FAMILY", envir=options()$htmlCounters)
  if (missing(size)) size <- get("HTML.FONT.SIZE", envir=options()$htmlCounters)
  
  
  if (reportFormat == "latex"){

    #*** Set title in bold text (can't be changed)    
    if (is.null(title)) title <- ""
    else title <- paste("\\textbf{", title, "}\\\\", sep="")

    #*** Set alignment command
    if (align %in% c("left", "right")) align <- paste("flush", align, sep="")
  
    size <- paste(map.size(size), "\n", sep="")

    #*** Opening String for alignment
    align.open <- paste("\\begin{", align, "}\n", sep="")

    #*** Closing String for alignment
    align.close <- paste("\\end{", align, "}", sep="")

    #*** Style String
    style <- ""
    if (italic)    style <- paste(style, "\\emph{", sep="")
    if (bold)      style <- paste(style, "\\textbf{", sep="")
    if (underline) style <- paste(style, "\\ul{", sep="")

    #*** Opening call to style
    style.open <- if (style %in% "") "%% \\emph{ \\textbf \\ul{\n" else style
    style.close <- paste(
                      paste(rep("}", italic + bold + underline), collapse=""),
                      sep="")
    style.close <- if (style.close %in% "") "%% } } } %%close emph, textbf, and ul" else style.close
    style.close <- paste(style.close, "\n", sep="")
           
    #*** Closing call to style            
    text <- paste(..., sep=sep)
    if (translate) text <- latexTranslate(text)
    text <- paste(text, "\n", sep="")
  
    #*** Paste all code together
    code <- paste(fncall, align.open, style.open, size, text,
                  style.close, align.close, "\n\n", sep="")
  }
  
  
  if (reportFormat == "html"){
    txt <- paste(..., sep=sep)
    
    code <- paste("<p style='font-family:", font, ", ", family, "; font-size:", map.size(size), "pt; ",
                  "text-align:", align, "; ",
                  "font-style:", if (italic) "italic; " else "none; ",
                  "font-weight:", if (bold) "bold; " else "none; ",
                  "text-decoration:", if (underline) "underline;" else "none;",
                  "'>", sep="")
    code <- paste(code, txt, "</p>")
  }
  
  if (reportFormat == "markdown"){
    code <- paste(..., sep=sep)
    code <- paste(if (bold) "**" else "", if (italic) "_" else "", code, if (italic) "_" else "", if (bold) "**" else "", sep="")
  }
  
  return(code)
}

