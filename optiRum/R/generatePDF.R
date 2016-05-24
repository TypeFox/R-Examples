#' Convert an .Rnw file to a PDF
#'
#' This function is designed to handle the production task of a 'standard' 
#' PDF process.  It is designed to build using pdflatex  (unless otherwise
#' specified) an adequate number of times to enable full typesetting to 
#' occur after taking into account items like contents pages, glossaries, 
#' and figures.
#'
#' @param srcpath     Location of .Rnw file, default is current directory
#' @param srcname     Rnw file name without extension e.g. 'Style'
#' @param destpath    Location of PDF file to be sent to, default is current directory
#' @param destname    PDF file name without extension e.g. 'Style_output'
#' @param DATED       Boolean indicating whether PDF filename should include yyyymmdd added to it
#' @param CLEANUP     Boolean indicating whether ancilliary files should be removed after production
#' @param QUIET       Boolean indicating whether console output should be limited
#' @param envir       Set default environment for knitr to run in - prevents a data.table issue
#' @param ...         Allows additional parameters to be passed to the knit2pdf function
#' 
#' @keywords knitr pdflatex generate PDF Rnw
#' @seealso \code{knit2pdf}
#' @family helper
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # simple call
#' generatePDF(srcname='basic')
#' 
#' # complex call
#' generatePDF(
#' srcname='basic',
#' destpath=getwd(),
#' destname='basic',
#' DATED=TRUE,
#' CLEANUP=FALSE,
#' QUIET=TRUE,
#' compiler='xelatex')
#' }
#' 
generatePDF <- function(srcpath = getwd(), srcname, destpath = getwd(), destname = srcname, DATED = FALSE, CLEANUP = TRUE, QUIET = FALSE, envir = new.env(parent = .GlobalEnv), 
    ...) {
    
    stopifnot(is.character(srcpath), is.character(srcname), is.character(destpath), is.character(destname), is.logical(DATED), is.logical(QUIET), is.logical(CLEANUP), 
        file.exists(file.path(srcpath, paste0(srcname, ".Rnw"))))
    
    src <- file.path(srcpath, paste0(srcname, ".Rnw"))
    dest <- file.path(destpath, paste0(destname, ifelse(DATED, format(Sys.Date(), "%Y%m%d"), ""), ".tex"))
    
    knitr::knit2pdf(input = src, output = dest, envir = envir, quiet = QUIET, clean = CLEANUP, ...)
    
} 
