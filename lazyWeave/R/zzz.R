.onLoad <- function(libname, pkgname){
  options("lazyReportFormat" = "markdown")
  options("lazyWeave_latexComments" = "markdown")
  options("lazyWeave_cat" = TRUE)
  options("htmlCounters" = new.env())
  assign("HTML.COUNTER.TABLE", 1, envir=options()$htmlCounters)
  assign("HTML.COUNTER.FIGURE", 1, envir=options()$htmlCounters)
  assign("HTML.COUNTER.FOOTNOTE", 1, envir=options()$htmlCounters)
  assign("HTML.COUNTER.CHAPTER", 1, envir=options()$htmlCounters)
  assign("HTML.COUNTER.SECTION", 1, envir=options()$htmlCounters)
  assign("HTML.COUNTER.SUBSECTION", 1, envir=options()$htmlCounters)
  assign("HTML.COUNTER.SUBSUBSECTION", 1, envir=options()$htmlCounters)
  assign("HTML.FOOTNOTES", NULL, envir=options()$htmlCounters)
  assign("HTML.FONT.FAMILY", "serif", envir=options()$htmlCounters)
  assign("HTML.FONT.FONT", "helvetica", envir=options()$htmlCounters)
  assign("HTML.FONT.SIZE", 11, envir=options()$htmlCounters)
}

