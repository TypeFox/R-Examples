#' @name lazy.file.start
#' @export lazy.file.start
#' 
#' @title Initiate New LaTeX, HTML, or Markdown Files
#' @description Write code to open a new LaTeX document with 
#'   packages, classes, a title, and page settings
#'   
#' @param docClass a character string giving a valid LaTeX document 
#' class. For example, \code{article}, \code{slide}, \code{report}, 
#' \code{book}
#' @param packages A character vector of additional LaTeX packages to use
#' @param counters A character vector of additional counters to initialize
#' @param layout LaTeX code for page layout.  Remember to escape backslashes!
#' @param page A character string denoting the page numbering style.  Options
#'   are \code{"arabic", "roman", "Roman", "alph", "Alph"}.
#' @param ligatures Determines if ligatures are enabled.  See the references for a link about ligatures
#' @param title A title for the document
#' @param author Author of the document
#' @param date Date to be printed on the title page
#' @param initialize For HTML and markdow files and when \code{TRUE}, 
#'   the function \code{lazy.options} is called and 
#'   all of the counters are reset to 1.  Font, family, and size 
#'   defaults are also reset 
#' 
#' @details Titles are only made when either \code{title} or \code{author} are 
#' not \code{NULL}.  
#' 
#' Packages automatically included are "xcolor", "graphicx", "colortbl", "float",
#' "soul", "hyperref", "placeins", and "Sweave".  Any user 
#' defined templates made in conjuction with \code{lazyWeave} must include these
#' packages in order to use figures and underlined text.
#' 
#' With \code{page}, the options produce the following:
#' \tabular{ll}{
#'   arabic \tab Arabic numbers\cr
#'   roman  \tab Lower case roman numerals (i, ii, iii, ...)\cr
#'   Roman  \tab Upper case roman numerals (I, II, III, ...)\cr
#'   alph   \tab Lower case alphabetic ordering (a, b, c, ...)\cr
#'   Alph   \tab Upper case alphabetic ordering (A, B, C, ...)\cr
#' }
#' 
#' @author Benjamin Nutter
#' 
#' @references Ligatures: \url{https://en.wikibooks.org/wiki/LaTeX/Text_Formatting#Ligatures}
#' 
#' @examples
#' 
#' #* lazy.file.start does not currently work with markdown documents
#' #* First, we set the lazyReportFormat option to "latex"
#' orig_option <- getOption("lazyReportFormat")
#' options(lazyReportFormat="latex")
#' lazy.file.start(docClass="report", 
#'   packages=c("pslatex", "palatino", "avant"),
#'   title="Report Name", author="Your Name")
#'  
#' #* Return the original option setting
#' options(lazyReportFormat=orig_option)
#'   

lazy.file.start <-
function(docClass="article", packages=NULL, 
    counters=NULL, layout="", page="arabic", ligatures=TRUE,
    title=NULL, author=NULL, date="", initialize=TRUE){
    
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html")) stop("option(\"lazyReportFormat\") must be either 'latex' or 'html'")
  
  
  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")
  
  
  if (reportFormat == "latex"){
    #*** Build a string for loading user designated packages    
    packages <- if (!is.null(packages)) paste(paste("\\usepackage{", packages, "}", sep=""), collapse="\n") else "%% \\usepackage{}"
      

    #*** Build a string to initialize counters
    counters <- if (!is.null(counters)) paste(paste("\\newcounter{", counters, "}", sep=""), collapse="\n") else "%% \\newcounter{}"
  
    #*** Build a string to set the title and author.  No title is made when both
    #*** title and author are NULL
    if (!is.null(title) | !is.null(author))
      title <- paste("\\title{", title, "}\n",
                     "\\author{", author, "}\n",
                     "\\date{", date, "}\n", 
                     "\\maketitle", sep="")
    else title <- ""

    #*** Paste all the elements together to open the LATEX file.
    code <- paste(fncall, "\n",
                  "\\documentclass{", docClass, "}\n",
                  "\\usepackage{breakurl, colortbl, fancyhdr, float, graphicx,\n", 
                  "             lastpage, lscape, microtype, soul, Sweave, url, xcolor}\n",
                  "\\usepackage[section]{placeins}\n",
                  packages, "\n",
                  if (!ligatures) "\\DisableLigatures{encoding=*, family=*}\n" else "",
                  if (layout %in% "") "%% layout commands may be written here" else layout, "\n",
                  counters, "\n",
                  "\\begin{document}\n\n",
                  lazy.page.number(page), "\n\n",
                  title, sep="")
  }
  
  #*** HTML Format
  #*** Yes, it does seem really lame, in comparison. But the HTML options
  #*** can be set using lazy.options.
  if (reportFormat == "html"){
    if (initialize){
      setHtmlOptions(table=1, figure=1, footnote=1, chapter=1, section=1, subsection=1, 
                     font.family="serif", font="helvetica", font.size=11)
      assign("HTML.FOOTNOTES", NULL, envir=options()$htmlCounters)
    }
    code <- "<html>\n"
  }
  
  if (reportFormat == "markdown"){
    setHtmlOptions(table=1, figure=1, footnote=1, chapter=1, section=1, subsection=1, 
                   font.family="serif", font="helvetica", font.size=11)
    assign("HTML.FOOTNOTES", NULL, envir=options()$htmlCounters)
    code <- ""
  }
  
  return(code)
}

