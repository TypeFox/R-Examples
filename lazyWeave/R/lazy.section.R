#' @name lazy.section
#' @export lazy.section
#' 
#' @title Open Sections in LaTeX Documents
#' @description Provides the code to start a new section, chapter, or
#'   subsection.
#' 
#' @param heading The name of the section or subsection
#' @param type Type of section.  \code{"section"} for section; 
#'   \code{"sub"} for subsection; and \code{"subsub"} for subsubsection.  
#'   The option \code{"sub2"} is a relic of me trying to avoid typing two 
#'   extra characters.  I decided it would be better to keep true to the
#'   LaTeX descriptors.  Thus, \code{"sub2"} is available for 
#'   back-compatibility, but I recommend against its use
#' @param ordered Logical.  Determines if the section is numbered
#' @param counter Name of a the counter to be used for this section.  
#'   When \code{NULL}, the value defaults to the counter corresponding
#'   to the type of section.  See \code{\link{lazy.counter}} for more 
#'   details about counters
#' @param counterSet Value to which \code{counter} should be set.  
#'   In other words, the number for this section (or similar division).
#' @param label The label to be used with \code{lazy.ref}.
#' @param font Font to be used in HTML files
#' @param family Font family to be used in HTML files
#' @param size Font size.  I'm pretty sure this doesn't actually get used, 
#'   but haven't gotten around to verifying this.  Heading sizes are set 
#'   using the HTML <H ... > tags
#' @param leadspace For HTML reports, should several lines of white space
#'   be placed before the section header.  This helps to create a visual
#'   break between sections
#' @param floatBarrier Figures and tables in LaTeX are floating objects 
#'   that LaTeX may choose to place somewhere other than where specified in
#'   order to optimize appearance.  This is not always desirable.  
#'   \code{floatBarrier} prevents floats from overlapping a section break.  
#'   This may be turned off if placement of figures and tables is of little 
#'   consequence
#' 
#' @details For HTML, Sections titles are printed using the \code{<Hx ...>} tags,
#'   where \code{x} is the heading number.  \code{H1} is the largest heading
#'   size and \code{H6} is the smallest.  Chapter headings use \code{<H2>}; 
#'   sections use \code{<H3>}; subsections use \code{<H4>}; and subsubsections 
#'   use \code{<H5>}.
#'   
#' @examples
#' \dontrun{
#'   lazy.write(
#'   lazy.file.start(),
#'   lazy.section("Section A", ordered=TRUE),
#'   lazy.text("Notice that Section A is numbered"),
#'   lazy.section("Subsection", type="sub", ordered=FALSE),
#'   lazy.text("But the subsection is not numbered"),
#'   lazy.file.end(),
#'   OutFile="Example-1.tex")
#'   
#'   unlink("Example-1.tex")
#'  }
#'  

lazy.section <- function(heading, type=c("section", "sub", "subsub", "chapter", "sub2"),
    ordered=FALSE, counter, counterSet=NULL, label=NULL,
    font=getOption("html.font.font"), family=getOption("html.font.family"),
    size=getOption("html.font.size"), leadspace=TRUE, floatBarrier=TRUE){
  
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex' 'html', or 'markdown'")
  
  #*** Construct the comment with the function call
  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")

  
  type <- match.arg(type, c("section", "sub", "subsub", "chapter", "sub2"))
        
  
  if (reportFormat == "latex"){
    counterStr <- if (!missing(counter)) lazy.counter(counter, fn="use") else "%% \\usecounter{}\n"
    if (!is.null(counterSet)) counterStr <- paste(counterStr, lazy.counter(counter, value=counterSet - 1, fn="set"), sep="\n")
  
    label <- if (!is.null(label)) lazy.label(label) else "%% \\label{}\n"

    if (type == "sub")                         type <- "subsection"
    if (type == "sub2" || type == "subsub")    type <- "subsubsection"
    if (type == "chapter")                     type <- "chapter"
  
    star <- if(ordered) "" else "*"
  
    code <- paste(fncall, counterStr, "\n", if (floatBarrier) "\\FloatBarrier\n" else "", 
                  "\\", type, star, "{", heading, "}\n", label, "\n\n", sep="")
  }
  
  
  if (reportFormat == "html"){
    if (missing(counter)) counter <- type
    
    if (!is.null(counterSet)) lazy.counter(counter, counterSet, fn="set")
    count.val <- lazy.counter(counter, fn="value")
    
    sec.number <- switch(type,
                         "chapter" = paste("Chapter ", count.val, ": ", sep=""),
                         "section" = paste(lazy.counter("chapter", fn="value") - 1, ".", count.val, ": ", sep=""),
                         "sub" = paste(lazy.counter("chapter", fn="value") - 1, ".", lazy.counter("section", fn="value") - 1, ".", count.val, ": ", sep=""),
                         "sub2" = paste(lazy.counter("chapter", fn="value") - 1, ".", lazy.counter("section", fn="value") - 1, ".",
                                        lazy.counter("sub", fn="value") - 1, ".", count.val, ": ", sep=""),
                         "subsub" = paste(lazy.counter("chapter", fn="value") - 1, ".", lazy.counter("section", fn="value") - 1, ".",
                                          lazy.counter("sub", fn="value") - 1, ".", count.val, ": ", sep=""))
    
    if (ordered){
      lazy.counter(counter, count.val + 1, fn="set")
      if (counter %in% "chapter"){
        lazy.counter("section", 1, "set")
        lazy.counter("sub", 1, "set")
        lazy.counter("subsub", 1, "set")
      }
      if (counter %in% "section"){
        lazy.counter("sub", 1, "set")
        lazy.counter("subsub", 1, "set")
      }
      if (counter %in% "sub") lazy.counter("subsub", 1, "set")
    }
    
    if (type %in% "section") H <- 3
    else if (type %in% "sub") H <- 4
    else if (type %in% c("subsub", "sub2")) H <- 5
    else H <- 2
    
    code <- paste(if (leadspace) "<br><br><br>" else "", #May need to change to </br> in the future.
                  "<H", H, " style='font-family:", font, ", ", family, "; font-weight:bold;'>",
                  if (ordered) sec.number else "",
                  heading, "</H", H, ">\n", sep="")  
  }
  
  if (reportFormat == "markdown"){
    if (missing(counter)) counter <- type
    
    if (!is.null(counterSet)) lazy.counter(counter, counterSet, fn="set")
    count.val <- lazy.counter(counter, fn="value")
    
    sec.number <- switch(type,
                         "chapter" = paste("Chapter ", count.val, ": ", sep=""),
                         "section" = paste(lazy.counter("chapter", fn="value") - 1, ".", count.val, ": ", sep=""),
                         "sub" = paste(lazy.counter("chapter", fn="value") - 1, ".", lazy.counter("section", fn="value") - 1, ".", count.val, ": ", sep=""),
                         "sub2" = paste(lazy.counter("chapter", fn="value") - 1, ".", lazy.counter("section", fn="value") - 1, ".",
                                        lazy.counter("sub", fn="value") - 1, ".", count.val, ": ", sep=""),
                         "subsub" = paste(lazy.counter("chapter", fn="value") - 1, ".", lazy.counter("section", fn="value") - 1, ".",
                                          lazy.counter("sub", fn="value") - 1, ".", count.val, ": ", sep=""))
    
    if (ordered){
      lazy.counter(counter, count.val + 1, fn="set")
      if (counter %in% "chapter"){
        lazy.counter("section", 1, "set")
        lazy.counter("sub", 1, "set")
        lazy.counter("subsub", 1, "set")
      }
      if (counter %in% "section"){
        lazy.counter("sub", 1, "set")
        lazy.counter("subsub", 1, "set")
      }
      if (counter %in% "sub") lazy.counter("subsub", 1, "set")
    }
    
    if (type %in% "section") H <- "##"
    else if (type %in% "sub") H <- "###"
    else if (type %in% c("subsub", "sub2")) H <- "####"
    else H <- "#"
    
    code <- paste(H, if (ordered) sec.number else "", heading, sep=" ")
    
    
  }
  
  return(code)
}
