#' @name lazy.list
#' @export lazy.list
#' 
#' @title Lists in LaTeX and HTML
#' @description Produce code for lists in LaTeX documents
#' 
#' @param item A vector with the items to be placed in the list
#' @param ordered Denotes if the list is ordered or bulleted
#' @param counter For future inclusion.  Specifies what counter 
#'   should be used for numbering.  Currently not in use
#' @param counterSet The value to which \code{counter} should 
#'   be set.  In other words, the number for the first item in the 
#'   list
#' @param title A title for the list
#' @param style A character string denoting how the ordered list 
#'   should be numbered.  Options are \code{"arabic", "roman", 
#'   "Roman", "alph", "Alph"}.
#' @param symbol A symbol for bulleted lists to be used as the bullet
#' @param font Font to be used in HTML documents.  Defaults to the 
#'   font set in the options
#' @param family Font family to be used in HTML documents.  Defaults 
#'   to the font family in the options
#' @param size Font size to be used in HTML documents.  Defaults to 
#'   the font family in the options
#'   
#' @details  With \code{style}, the options produce the following and 
#' apply to both LaTeX and HTML:
#' \tabular{ll}{
#'   arabic \tab Arabic numbers\cr
#'   roman  \tab Lower case roman numerals (i, ii, iii, ...)\cr
#'   Roman  \tab Upper case roman numerals (I, II, III, ...)\cr
#'   alph   \tab Lower case alphabetic ordering (a, b, c, ...)\cr
#'   Alph   \tab Upper case alphabetic ordering (A, B, C, ...)\cr
#' }
#' 
#' The given options for \code{symbol} follow the HTML conventions, 
#' but when \code{options("lazyReportFormat")} is \code{"latex"}, 
#' these options are translated into the appropriate equivalent.
#' 
#' @author Benjamin Nutter
#' 
#' @examples
#' \dontrun{
#' lazy.write(
#'   lazy.file.start(),
#'   lazy.text("A vector can be presented as a list as follows"),
#'   lazy.list(c("First item", "Second item", "Third item", 
#'               "Fourth item", "Fifth item"), style="Alph"),
#'   lazy.file.end(),
#'   OutFile="Example 1.tex")
#' 
#' unlink("Example 1.tex")
#' }
#' 

lazy.list <-
function(item, ordered=TRUE, counter=NULL, counterSet=1, title=NULL, 
         style=c("arabic", "Roman", "roman", "Alph", "alph"), 
         symbol=c("bullet", "circ", "blacksquare"),
         font, family, size){
  #$\\bullet$, $\\circ$, $\\blacksquare$,
  
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  style.ref <- data.frame(latex=c("arabic", "Roman", "roman", "Alph", "alph"),
                          html =c("arabic", "I",     "i",     "A",    "a"),
                          stringsAsFactors=FALSE)
  symbol.ref <- data.frame(latex=c("bullet", "circ", "blacksquare"),
                           html =c("circle", "disc", "square"),
                           stringsAsFactors=FALSE)
  
  if (missing(font)) font <- get("HTML.FONT.FONT", envir=options()$htmlCounters)
  if (missing(family)) family <- get("HTML.FONT.FAMILY", envir=options()$htmlCounters)
  if (missing(size)) size <- get("HTML.FONT.SIZE", envir=options()$htmlCounters)
  
  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")  
  
  #*** Options for style
  #*** arabic, Roman, roman, Alph, alph
  style <- match.arg(style, c("arabic", "Roman", "roman", "Alph", "alph"))
  symbol <- match.arg(symbol, c("bullet", "circ", "blacksquare"))
  
  if (reportFormat == "latex"){
    symbol <- paste("$\\", symbol, "$", sep="")

    #*** Print the title of the list in bold face (currently there is no way 
    #*** to turn off boldface type
    if (is.null(title)) title <- "%% \\textbf{} %% List title" 
    else{ 
      title <- paste("\\textbf{", title, "}", sep="")
    }

    #*** Make a string of the items in the list
    item <- paste(paste("\\item", item), collapse="\n  ")
    item[1] <- paste(" ", item[1])

    #*** If no counter for an ordered list is given, this makes a new counter
    if (is.null(counter)){
      newcount <- paste(sample(LETTERS, 6), collapse="")
      code <- lazy.counter(newcount)
    }
    else code <- "%% \\newcounter{}\n"
  
    #*** make string for an ordered list, else unordered list
    if (ordered){
      orderedStart <- lazy.counter(if (is.null(counter)) newcount else counter, value=counterSet - 1, fn="set")
      if (is.null(counter)){
        lst <- paste("\\begin{list}{\\", style, "{", newcount, "}}",
                     "{\\usecounter{", newcount, "}}", sep="")
      }
      else{
        lst <- paste("\\begin{list}{\\", style, "{", counter, "}}",
                     "{\\usecounter{", counter, "}}", sep="")
      } 
      code <- paste(code, title, lst, orderedStart, sep="\n")
    }
    else{
      orderedStart <- ""
      code <-paste(title, "\n\\begin{list}{", symbol, "}{}", sep="")
    }

    #*** Paste code together for the list
    code <- paste(fncall, code, item, "\\end{list}", "\n\n", sep="")
  }
  
  if (reportFormat == "html"){
    #*** Match style and symbol from latex arguments to html arguments
    style <- style.ref[style.ref$latex == style, "html"]
    symbol <- symbol.ref[symbol.ref$latex == symbol, "html"]
    if (length(symbol) == 0) symbol <- "circle"
    
    #if (is.null(counter)) lazy.counter("html.counter.list", fn="set", value=1)

    #*** Print the title of the list in bold face (currently there is no way 
    #*** to turn off boldface type
    if (is.null(title)) title <- "<!-- List title -->" 
    else{ 
      title <- lazy.text.format(title, bold=TRUE)
    }
    
    tag <- if (ordered) "ol" else "ul"
    type <- if (ordered) style else symbol
    if (type %in% "arabic") type <- ""
    
    code <- paste("<", tag, 
                  " start='", if (is.null(counter)) 1 else lazy.counter(counter, fn="value"), 
                  "' type='", type, "' ",
                  "style='font-family:", font, ", ", family, "; font-size:", size, "pt;'>", title, sep="")
    lst <- paste("  <li>", item)
    lst <- paste(lst, collapse="\n")
    
    code <- paste(fncall, code, "\n", lst, "\n</", tag, ">\n\n", sep="")  
  }
  
  if (reportFormat == "markdown"){
    if (!ordered) code <- paste(paste("*", item), collapse="\n")
    else {
      val <- 1:length(item) + (counterSet-1)
      code <- paste(paste(val, ". ", item, sep=""), collapse="\n")
    }
  }

  return(code)
}

