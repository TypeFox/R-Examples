#' @name lazy.write
#' @export lazy.write
#' 
#' @title Output LaTeX Code to a File
#' @description Output text and LaTeX code to a file for building
#' 
#' @param ... Strings, expressions, and statements to be written 
#'   to a .tex file
#' @param OutFile Filename to which the code should be written
#' @param append logical.  Indicates if the code should be appended to 
#'   and existing file
#' @param collapse Sets the character string that is placed 
#'   between elements in \code{\dots}
#' @param footnotes logical.  For HTML and markdown only, when 
#'   \code{TRUE}, the footnotes stored in 
#'   \code{options("lazy.footnotes")} will be appended to the end 
#'   of the HTML document
#'   
#' @details The contents of \code{\dots} will be pasted together
#' 
#' @author Benjamin Nutter
#' 

lazy.write <-
function(..., OutFile, append=FALSE, collapse="\n", footnotes=TRUE){
  
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")

  if (reportFormat != "markdown"){
  #*** Stop Function if no OutFile is given
    if (missing(OutFile))
      stop("'OutFile' must be explicitly specfied using 'OutFile=[filename]'")
  

      file <- unlist(strsplit(OutFile, "[.]"))
      file.ext <- utils::tail(file, 1)
      if (reportFormat == "latex" && file.ext %in% c("html", "htm")) 
        OutFile <- paste(file[-length(file)], ".tex", sep="")
      if (reportFormat == "html" && file.ext == "tex")
        OutFile <- paste(file[-length(file)], "html", sep=".")
  }

  #*** As awful as this sounds, I have no idea what this does.  It is a relic
  #*** from html.write in my earlier version of CCFmisc, but I never documented
  #*** why I needed this.  I think it has something to do with making vectors
  #*** appear across a row instead of as a column.  But until the function 
  #*** shows behavior contrary to what I desire, I don't want to bother 
  #*** investigating what this is doing.
  f <- function(x){
    x <- ifelse(!is.null(dim(x)) || is.vector(x), paste(t(x), collapse=" "), x)
    return(x)
  }


  #*** Combine all the code into one string for exporting to the file.
  code <- list(...)
  code <- if (!is.null(code)) lapply(code, f) else code <- ""
  code <- paste(code, collapse = collapse)
  
  if (reportFormat == "markdown" & footnotes & !is.null(get("HTML.FOOTNOTES", envir=options()$htmlCounters))){
    return(get("HTML.FOOTNOTES", envir=options()$htmlCounters))
  }
  else if (reportFormat=="markdown") return("")
  
  if (reportFormat == "html" & footnotes & !is.null(get("HTML.FOOTNOTES", envir=options()$htmlCounters))){
    code <- gsub("</html>", paste("\n\n\n", get("HTML.FOOTNOTES", envir=options()$htmlCounters), "\n</html>"), code)
    assign("HTML.FOOTNOTES", NULL, envir=options()$htmlCounters)
  }
  
  write(code, OutFile, append = append)
}

