#' @name lazy.figure
#' @export lazy.figure
#' 
#' @title Include Figures in Latex Documents
#' @details Generates the code to place a figure in a Latex document
#' 
#' @param filename Character string giving the location of the file to be included
#' @param caption Text giving the caption for the figure
#' @param align Character string stating the alignment.  Valid options are
#'   \code{"left"}, \code{"right"}, or \code{"center"}
#' @param height The height of the figure
#' @param width The width of the figure
#' @param units The units for height and width.  Defaults to \code{"in"}.
#' @param counter Name of a counter to use to number the table
#' @param counterSet The number to which \code{counter} should be set.  
#'   In other words, the figure number for this figure
#' @param label Name of a label
#' @param placement Controls the placement of the figure.  Options are
#'   \code{"ht", "t", "b", "p", "H"} and can be supplemented with 
#'   \code{"!"}. See "Details" for more explanation
#' @param alt For HTML documents only--when \code{filename} cannot be found, 
#'   this text is printed in the figure's place
#' @param cat Logical. Determines if the output is returned as a character string
#'   or returned via the \code{cat} function (printed to console).  The default
#'   value is set by \code{options()$lazyWeave_cat}.  This argument allows for
#'   selective override of the default.
#' 
#' @details
#' For LaTeX files, \code{placement} options are used as follows:
#' \tabular{ll}{
#' ht \tab Place the float here, i.e., 
#' approximately at the same point it occurs \cr
#' t  \tab Position at the top of the page\cr
#' b  \tab Position at the bottom of the page \cr
#' p  \tab Put on a special page for floats only \cr
#' H  \tab Places the float at precisely the location in the LaTeX code. 
#' Requires the float package\cr
#' }
#' The \code{"!"} may be used after any of these in order to override 
#' LaTeX float rules and force your selection.  More can be learned by 
#' reading about floats in a LaTeX manual.
#' 
#' For HTML files, the file can be any type supported by the browser.  JPEGs and PNGs seem to work well.
#' 
#' @author Benjamin Nutter
#' @examples
#' \dontrun{
#' pdf("MPG.pdf", height=4, width=4)
#' hist(mtcars$mpg)
#' dev.off()
#' 
#' lazy.figure("MPG.pdf")
#' 
#' lazy.write(
#' lazy.file.start(),
#' lazy.figure("MPG.pdf", 
#' caption="Distribution of Miles per Gallon in mtcars dataset",
#' height=5, width=5, label="MPGgraph"),
#' lazy.file.end(),
#' OutFile="Example-1.tex")
#' 
#' unlink("MPG.pdf")
#' unlink("Example-1.tex")
#' unlink("Example-1.pdf")
#' }

lazy.figure <-
function(filename, caption=NULL, align="center",
                         height=3, width=3, units="in", 
                         counter, counterSet=NULL,
                         label=NULL, placement="h",
                         alt="Image Not Found",
         cat=getOption("lazyWeave_cat")){

  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  
  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")
  
  #*** LaTeX format
  if (reportFormat == "latex"){

    #*** Set align argument to latex command
    if (align %in% c("left", "right")) align <- paste("flush", align, sep="")

    #*** Set height and width strings for figure
    height <- paste("height=", height, units, sep="")
    width  <- paste("width=",  width, units, sep="")
   
    #*** Specify and set counter
    counterStr <- if (!missing(counter)) paste("  ", lazy.counter("figure", fn="use"), sep="") else "%% \\usecounter{}\n"
    if (!is.null(counterSet)) counterStr <- paste(counterStr, "\n  ", lazy.counter(counter, value=counterSet - 1, fn="set"), sep="")

    #*** Set caption and label strings
    caption <- if (is.null(caption))  "      %% \\caption{}\n"
      else paste("      \\caption{", caption, "}\n", sep="")
    label <- if (is.null(label))  "      %% \\label{}\n"
      else paste("      \\label{", label, "}\n", sep="")

    #*** Produce LATEX code for the figure.
    code <- paste(
      fncall,
      "\\begin{figure}[", placement, "]\n",
      counterStr,
      "  \\begin{", align, "}\n",
      "    \\includegraphics[",height, ", ", width, "]{", filename, "}\n",
      caption, 
      label, 
      "  \\end{", align, "}\n",
      "\\end{figure}", sep="")
  }
  
  
  #*** HTML format
  if (reportFormat == "html"){
    
    if (missing(counter)) counter <- "figure"
    #*** Caption
    if (is.null(caption)) caption <- ""
    else{
      if (!is.null(counterSet)) lazy.counter(counter, counterSet, fn="set")
      count.val <- lazy.counter(counter, fn="value")
      caption <- paste("Figure ", lazy.counter(counter, fn="value"), ": ", caption, sep="")
      lazy.counter(counter, count.val + 1, fn="set")
    }
    
    code <- paste(fncall,
                  "<p style='text-align:", align, ";'>",
                  "<img src='", filename, "' height=", height, units, " width=", width, units, " alt='", alt, "'/></p>\n", 
                  lazy.text(caption, italic=TRUE, align=align), sep="")
    if (!is.null(label)) code <- paste(lazy.label(label), code, sep="\n") 
    code <- paste(code, "\n\n", sep="")
  }
  
  #*** Markdown format
  if (reportFormat == "markdown"){
  
    if (missing(counter)) counter <- "figure"
    if (is.null(caption)) caption <- ""
    else{
      if (!is.null(counterSet)) lazy.counter(counter, counterSet, fn="set")
      count.val <- lazy.counter(counter, fn="value")
      caption <- paste("Figure ", lazy.counter(counter, fn="value"), ": ", caption, sep="")
      lazy.counter(counter, count.val + 1, fn="set")
    }
    
    code <- paste("![", alt, "][", lazy.counter(counter, fn="value"), "]", "\n\n",
                  "[", lazy.counter(counter, fn="value"), "]: ", filename, " ''", sep="")
    
    
  }

  if (cat) cat(code)
  else return(code)
}
