#' Generate Latex or HTML reports from R
#'
#' Arguments passed to functions are embedded in Latex code and can
#' be output to a file.  Allows the use of Latex and HTML to write reports, functions
#' for building tables, etc.  
#'
#' Depending on the working directory (for Windows users), the user may encounter the error 
#' "texi2dvi.exe: Windows API error 5: Access is denied." when trying
#' to build documents.  This happens when the working directory is write protected.  It is advisable to change the working directory to 
#' something not write protected.  This can be done permanently by right clicking on the R shortcut icon, selecting properties, and changing
#' the directory in the "Start in:" box.
#' 
#' \code{lazyWeave} assumes the availability of the packages \code{xcolor, graphicx, colortbl, soul, lscape,} and \code{Sweave}.  
#' If these packages are not available, the package most likely will not function properly. 
#' 
#' It should be noted that \code{lazyWeave} is a rather inefficient way to go 
#' about writing reports with LaTeX or HTML.  It's only real advantage is it reduces the
#' amount of knowledge a user needs to have about LaTeX (and it could be 
#' debated if that's really an advantage). 
#' 
#' Use of \code{lazyWeave} could also be greatly supplemented by some basic
#' familiarity with LaTeX.  For example, knowing the commands for bolding 
#' (\\textbf\{\}), italicizing (\\emph\{\}), and underlining text 
#' (\\ul\{\}) can go a long way to improving the look of reports.  It also would
#' help to know how to subscript and superscript terms.  Most introductions
#' to LaTeX will cover these basics. 
#' 
#' \code{lazyWeave} is also only intended to provide the most basic functionality
#' of LaTeX, and I have no plans of extending it much further than what is 
#' available already.  If what is in the package now is not sufficient enough 
#' to satisfy your needs, then I strongly suggest you look into using 
#' \code{Sweave}.
#' 
#' All of the functions can be used for LaTeX and HTML reports, but the functionality and appearance
#' may not be identical between formats.  
#'
#' @name lazyWeave
#' @docType package

NULL
