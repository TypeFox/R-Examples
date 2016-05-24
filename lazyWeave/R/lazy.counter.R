#' @name lazy.counter
#' @export lazy.counter
#' 
#' @title Create and Manage Counters for LaTeX Documents
#' @description Provides the code to create, manipulate, or extract values 
#'   of counters in a LaTeX or HTML document
#'   
#' @param counter A character(1) giving the name of the counter to be created
#'   and/or manipulated
#' @param value An integer.  For \code{fn="addto"}, it is the value by which 
#'   the counter is to be increased.  For \code{fn="set"}, it is the value to 
#'   which the counter should be set
#' @param oldcounter character(1).  An optional argument for \code{fn="new"}.  
#'   It must be a previously named counter.  If present, the new counter will 
#'   be reset whenever \code{oldcounter} is incremented 
#' @param fn Selects the LaTeX function to be used
#' 
#' @details
#' Counters are used to provide table, figure, and section numbers.  After each use of each command, the counter is incremented so that it
#' can be referred to in later uses.  New counters may be defined by users, but several LaTeX environments have default counters that do
#' not need to be defined.  These are \code{part, chapter, section, subsection, subsubsection, paragraph, subparagraph, page, equation,
#'                                          figure, table, footnote, mpfootnote}.  Any of these may be manipulated by \code{lazyWeave}.
#' 
#' Referring to and manipulating counters is done using \code{lazy.counter}.  Different actions are achieved by changing the \code{fn}
#' argument.  
#' 
#' \code{fn="new"} creates a new counter with name \code{counter}.
#' 
#' \code{fn="addto"} adds \code{value} to the current value of \code{counter}.
#' 
#' \code{fn="set"} changes the current value of \code{counter} to \code{value}.
#' 
#' \code{fn="use"} designates \code{counter} for use in the current environment.
#' 
#' \code{fn="value"} returns the value of \code{counter}.  This value isn't printed, but can be useful for doing arithmetic with 
#' counter values.
#' 
#' The HTML counters in \code{counter} are the defaults and will reference global counters 
#' \tabular{l}{
#'  \code{options("html.counter.table")}\cr
#'  \code{options("html.counter.table"}\cr
#'  \code{options("html.counter.table"}\cr
#'  \code{options("html.counter.table"}\cr
#'  \code{options("html.counter.table"}\cr
#'     \code{options("html.counter.table"}\cr
#'     \code{options("html.counter.table"}
#'   }
#'   
#'   Additional and custom counters may be defined if desired, in which case a new option will be defined as
#'   \code{options("html.custom.[countername]")}.
#'   
#'   Extracting a current value does not increment the counter--this must be done manually (and is done manually when used in \code{lazy.table}, \code{lazy.figure},\
#'   \code{lazy.footnote}, and \code{lazy.section}.
#'
#' @author Benjamin Nutter
#' 
#' @examples
#' lazy.counter("myCounter")
#' lazy.counter("myCounter", value=3, fn="set")
#' lazy.counter("myCounter", value=2, fn="addto")
#' lazy.counter("myCounter", fn="use")
#' lazy.counter("myCounter", fn="value")
#' 
#' lazy.counter("table", fn="use")

lazy.counter <- function(counter, value, oldcounter, fn=c("new", "addto", "set", "use", "value")){

  fn <- match.arg(fn, c("new", "addto", "set", "use", "value"))

  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")

  
  #*** Latex Counters
  if (reportFormat == "latex"){
    if (!missing(value)) if (!is.numeric(value)) stop("'value' must be numeric")
  
    #*** newcounter function
    if (fn %in% "new"){
      txt <- paste("\\newcounter{", counter, "}", sep="")
      if (!missing(oldcounter)) txt <- paste(txt, "[", oldcounter, "]", sep="")
    }
  
    #*** addtocounter function
    if (fn %in% "addto"){
      txt <- paste("\\addtocounter{", counter, "}{", value, "}", sep="")
    }
  
    #*** setcounter function
    if (fn %in% "set"){
      txt <- paste("\\setcounter{", counter, "}{", value, "}", sep="")
    }
  
    #*** usecounter function
    if (fn %in% "use"){
      txt <- paste("\\usecounter{", counter, "}", sep="")
    }
  
    #*** value function
    if (fn %in% "value"){
      txt <- paste("\\value{", counter, "}", sep="")
    }
  
    return(txt)
  }
  
  
  #*** HTML counters
  if (reportFormat %in% c("html", "markdown")){
    if ("set" %in% fn){
      if ("table" %in% counter)         assign("HTML.COUNTER.TABLE", value, envir=options()$htmlCounters)
      else if ("figure" %in% counter)   assign("HTML.COUNTER.FIGURE", value, envir=options()$htmlCounters)
      else if ("footnote" %in% counter) assign("HTML.COUNTER.FOOTNOTE", value, envir=options()$htmlCounters)
      else if ("chapter" %in% counter)  assign("HTML.COUNTER.CHAPTER", value, envir=options()$htmlCounters)
      else if ("section" %in% counter)  assign("HTML.COUNTER.SECTION", value, envir=options()$htmlCounters)
      else if ("sub" %in% counter)      assign("HTML.COUNTER.SUBSECTION", value, envir=options()$htmlCounters)
      else if ("sub2" %in% counter || "subsub" %in% counter)    assign("HTML.COUNTER.SUBSUBSECTION", value, envir=options()$htmlCounters)
      else assign(paste("HTML.COUNTER.", counter, sep=""), value, envir=options()$htmlCounters)
    }
    else if ("new" %in% fn){
      assign(paste("HTML.COUNTER", counter, ".", sep=""), if (missing(value)) 1 else value, envir=options()$htmlCounters)
    }
    else if ("value" %in% fn){
      if (counter %in% "table")    return(get("HTML.COUNTER.TABLE", envir=options()$htmlCounters))
      else if (counter %in% "figure")   return(get("HTML.COUNTER.FIGURE", envir=options()$htmlCounters))
      else if (counter %in% "footnote") return(get("HTML.COUNTER.FOOTNOTE", envir=options()$htmlCounters))
      else if (counter %in% "chapter")  return(get("HTML.COUNTER.CHAPTER", envir=options()$htmlCounters))
      else if (counter %in% "section")  return(get("HTML.COUNTER.SECTION", envir=options()$htmlCounters))
      else if (counter %in% "sub")      return(get("HTML.COUNTER.SUBSECTION", envir=options()$htmlCounters))
      else if (counter %in% "sub2")     return(get("HTML.COUNTER.SUBSUBSECTION", envir=options()$htmlCounters))
      else get(paste("HTML.COUNTER.", counter, sep=""), envir=options()$htmlCounters)

    }
    else message("The functions 'addto' and 'use' are not defined for HTML format.  No action is taken.")
  }
}
  
  
